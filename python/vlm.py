"""Zero-shot event classification through a vision-language model.

Each event's per-strip trace is rasterized to a small RGB image and pushed
through a VLM, which is asked which of four classes it is: (a,n), (a,a'),
beam, or none of those. Two things make this cheap enough to run over a
million events:

- **Nothing touches the disk and nothing calls matplotlib.** render_traces
  builds the whole batch of images with numpy array ops (~microseconds an
  event); the arrays are handed to the processor as in-memory PIL images.
  A matplotlib figure costs 10-30 ms, which over the blind reservoir would
  be hours of CPU spent before the GPU sees anything.

- **The model never generates.** class_probs runs ONE forward pass and reads
  the logits at the last prompt position, softmaxed over the four class token
  ids alone. That is the same cost as emitting a single token, it removes any
  dependence on the model obeying a formatting instruction (so a small model
  that would otherwise ramble still scores cleanly), and it yields a
  posterior per class rather than a hard label -- which is what the
  downstream unfolding wants, since a threshold can then be swept instead of
  taking whatever the model asserted.

Prefill cost is linear in the visual token budget, so config.VLM_TOKEN_BUDGET
is the knob that sets the run time. The budget the processor ACTUALLY
produced is measured at construction rather than assumed.

The processor half lives in PromptBuilder and the model half in
EventClassifier, deliberately: everything that can go wrong with the chat
template, the padding side, the visual token budget and the class token ids
is then checkable without downloading 8B of weights. vlm_selftest.py does
exactly that before a full run is worth starting.

Run inside the dataset dev shell, from python/. Gemma is a gated repo on the
Hub: accept the license and have a token in the environment (HF_TOKEN) or in
~/.cache/huggingface before the first run.
"""

import os
import time

import numpy as np

import config

# Set before torch is imported (it is imported lazily inside the classes, so
# this wins). Gemma 4's vision tower allocates one very large short-lived
# one_hot per forward; with the default allocator those come back as
# unusable reserved-but-unallocated blocks and a batch that fits in
# principle fails in practice. Overridable from the environment.
os.environ.setdefault("PYTORCH_CUDA_ALLOC_CONF", "expandable_segments:True")

# Candidate names for the Gemma 4 image processor's soft-token budget. The
# model card documents the budget but not the argument that sets it, so every
# plausible spelling is tried and the result is MEASURED (see
# PromptBuilder.measure_budget) rather than trusted.
#
# Where the Gemma 4 image processor keeps its soft-token budget: confirmed by
# the stage-2 dump to be Gemma4ImageProcessor.max_soft_tokens, default 280.
# It is a CEILING, not a count -- the processor picks the largest patch grid
# fitting under it, so 280 yields 256 (16x16) because 17x17 = 289 would not
# fit, and small source images are resized up to use the allowance.
#
# It applies to EVERY image in the call, so the reference figure and the
# event share one allowance and cannot be budgeted apart in a single
# processor call. The reference is byte-identical on every request and leads
# the prompt, so its share is the same tokens recomputed once per event --
# the strongest argument for caching its KV and expanding across the batch,
# if this is ever worth optimizing.
# Inputs indexed by IMAGE rather than by token. Both must be sliced together
# when the prompt is split into a cached prefix and a per-event tail.
_IMAGE_INDEXED = ("pixel_values", "image_position_ids")

_BUDGET_ATTRS = ("max_soft_tokens", "token_budget", "visual_token_budget",
                 "num_image_tokens", "image_token_budget", "max_image_tokens")


AN_INDEX = None  # resolved on first use; index of the (a,n) column


def an_probability(probs):
    """The (a,n) posterior column of a (n, n_classes) probability array.

    The only column any decision rests on: (a,a') and beam are not results,
    just somewhere for non-(a,n) events to go instead of contaminating the
    one column that feeds a cross section."""
    global AN_INDEX
    if AN_INDEX is None:
        AN_INDEX = list(config.VLM_CLASSES).index("an")
    return probs[:, AN_INDEX]


def is_an(probs, threshold=None):
    """Boolean (a,n) tag at a threshold on p(an), default
    config.VLM_AN_THRESHOLD. Deliberately not the argmax -- see that knob."""
    thr = config.VLM_AN_THRESHOLD if threshold is None else threshold
    return an_probability(probs) >= thr


def seed_perturbation(seed):
    """The (letter_order, dy, dhalf) this seed prompts and renders with.

    Perturbing these is what makes a seed ensemble mean anything: the forward
    pass is deterministic and the label is an argmax, so a seed that changed
    only a subsample would leave the answers bit-identical and report a
    systematic of exactly zero. Both axes are choices the method makes with
    no physics behind them -- which letter stands for which class, and where
    the polyline lands to the nearest fraction of a pixel.

    The first seed in config.VLM_SEEDS is deliberately unperturbed, so one
    run in the ensemble is always the canonical image and prompt.
    """
    n_cls = len(config.VLM_CLASSES)
    order = tuple(range(n_cls))
    dy = dhalf = 0.0
    if seed == config.VLM_SEEDS[0]:
        return order, dy, dhalf
    rng = np.random.default_rng(seed)
    if config.VLM_SEED_PERMUTE_LABELS:
        order = tuple(int(i) for i in rng.permutation(n_cls))
    if config.VLM_SEED_JITTER:
        dy = float(rng.uniform(-config.VLM_JITTER_DY, config.VLM_JITTER_DY))
        dhalf = float(
            rng.uniform(-config.VLM_JITTER_DHALF, config.VLM_JITTER_DHALF))
    return order, dy, dhalf


def build_prompt(letter_order, reference_strip=None):
    """The prompt, with class descriptions assigned to letters.

    `letter_order` is a permutation of indices into config.VLM_CLASSES:
    position i gets config.VLM_CLASS_TOKENS[i], so letter_order[i] names the
    class that letter stands for. Composed rather than hard-coded because
    which letter means which class is one of the arbitrary choices the seed
    ensemble perturbs -- VLMs carry label-token and option-order biases, so
    the assignment is a systematic, not a formatting detail.

    `reference_strip`, when given, adds the paragraph that explains the
    worked-example figure shown ahead of the event and says which reaction
    strip it was drawn for. Naming the strip is deliberate: the whole point
    of the test is whether one strip's example generalizes to the others.
    """
    lines = [config.VLM_PROMPT_HEADER]
    if reference_strip is not None:
        lines.append(config.VLM_REFERENCE_DESC.format(strip=reference_strip))
    lines.append(config.VLM_PROMPT_QUESTION)
    for i, cls in enumerate(letter_order):
        name = config.VLM_CLASSES[cls]
        lines.append(f"{config.VLM_CLASS_TOKENS[i]}: "
                     f"{config.VLM_CLASS_DESC[name]}")
    lines.append(config.VLM_PROMPT_FOOTER)
    return "\n".join(lines)


def _load_reference():
    """The worked-example figure as a PIL image, or None.

    config.VLM_REFERENCE_IMAGE is a plot the C++ already writes, so it is
    absent until strip-sum-scatter has run. That is not fatal -- the prompt
    falls back to its text descriptions -- but it is announced, because
    silently dropping the example would change what is being measured.
    """
    path = config.VLM_REFERENCE_IMAGE
    if path is None:
        return None
    if not path.is_file():
        print(f"vlm: no reference figure at {path}; prompting from the text "
              "descriptions alone")
        return None
    from PIL import Image
    return Image.open(path).convert("RGB")


def _hf_token():
    """Hub token for the gated Gemma repos, or None.

    $HF_TOKEN wins when it is set; otherwise config.VLM_HF_TOKEN_FILE is read
    if it exists. None is returned rather than raised, so a cached
    `hf auth login` still works and the Hub produces its own error if the
    repo really is unreachable."""
    env = os.environ.get("HF_TOKEN")
    if env:
        return env.strip()
    path = config.VLM_HF_TOKEN_FILE
    if path.is_file():
        return path.read_text().strip() or None
    return None


def _polyline_mask(values, height, width, lo, hi, half, dy=0.0):
    """Boolean (n, height, width) mask of the polyline through `values`.

    `values` is (n, m) -- m samples per row, one per strip. Columns are the
    linear interpolation of those samples across `width`; rows run top-down
    with `hi` at row 0, so the image reads like a plot. Adjacent columns are
    connected by filling the span between neighbouring samples, which is what
    keeps a steep step from breaking into disconnected dots, and `half`
    thickens the line by that many rows either side.

    Values outside [lo, hi] are clipped just off the canvas, so an off-scale
    excursion leaves the frame instead of piling up on the edge row.
    """
    m = values.shape[1]
    u = np.linspace(0.0, m - 1.0, width)
    i0 = np.floor(u).astype(np.int64)
    i1 = np.minimum(i0 + 1, m - 1)
    frac = (u - i0)[np.newaxis, :]
    y = values[:, i0] * (1.0 - frac) + values[:, i1] * frac
    row = (hi - y) / (hi - lo) * (height - 1.0) + dy
    row = np.clip(row, -(half + 1.0), height + half)
    prev = np.concatenate([row[:, :1], row[:, :-1]], axis=1)
    nxt = np.concatenate([row[:, 1:], row[:, -1:]], axis=1)
    top = np.minimum(row, np.minimum(prev, nxt)) - half
    bot = np.maximum(row, np.maximum(prev, nxt)) + half
    rows = np.arange(height, dtype=np.float64)[np.newaxis, :, np.newaxis]
    return (rows >= top[:, np.newaxis, :]) & (rows <= bot[:, np.newaxis, :])


def render_traces(X, beam_ref=None, X2=None, dy=0.0, dhalf=0.0):
    """Rasterize per-strip traces to (n, H, W, 3) uint8 RGB images.

    `X` is (n, m) in calibrated a.u. (beam ~ 1). The y-axis is FIXED at
    config.VLM_DE_MIN..VLM_DE_MAX for every event -- absolute, never
    per-trace max-normalized, because the amplitude above beam is part of
    what separates (a,a') from (a,n).

    `beam_ref` (m,) draws the unreacted-beam reference underneath the trace,
    which gives the model an in-image baseline instead of asking it to
    remember one. `X2` (n, m) draws a second curve in a third colour, for the
    L/R split view (config.VLM_SPLIT_LR) where the two ends carry information
    their sum does not.

    `dy` shifts every curve by that many rows and `dhalf` widens or
    narrows the line: the seed ensemble's rasterization jitter. Both default
    to no perturbation, so an unseeded call renders the canonical image.

    Pure numpy over the whole batch: no figure, no file, no encode.
    """
    X = np.asarray(X, dtype=np.float64)
    n, m = X.shape
    h, w = config.VLM_IMG_H, config.VLM_IMG_W
    lo, hi = config.VLM_DE_MIN, config.VLM_DE_MAX
    half = max(config.VLM_LINE_HALFWIDTH + dhalf, 0.0)
    img = np.empty((n, h, w, 3), dtype=np.uint8)
    img[:] = np.asarray(config.VLM_BG_RGB, dtype=np.uint8)
    # Reference first, then the second curve, then the trace on top, so the
    # event itself is never occluded by what it is being compared against.
    if beam_ref is not None and config.VLM_DRAW_BEAM_REF:
        ref = np.broadcast_to(
            np.asarray(beam_ref, dtype=np.float64)[:m], (n, m))
        img[_polyline_mask(ref, h, w, lo, hi, half,
                           dy)] = config.VLM_BEAM_RGB
    if X2 is not None:
        second = np.asarray(X2, dtype=np.float64)
        img[_polyline_mask(second, h, w, lo, hi, half,
                           dy)] = config.VLM_SECOND_RGB
    img[_polyline_mask(X, h, w, lo, hi, half, dy)] = config.VLM_TRACE_RGB
    return img


class PromptBuilder:
    """The processor half: chat template, token budget, class token ids.

    Holds no model, so it can be constructed from the Hub's processor files
    alone -- a few MB rather than the full checkpoint. That is what makes the
    templating independently testable (vlm_selftest.py stage 2) instead of
    only failing an hour into a run.
    """

    def __init__(self, model_id=None, token_budget=None, letter_order=None):
        from transformers import AutoProcessor

        self.model_id = model_id or config.VLM_MODEL
        n_cls = len(config.VLM_CLASSES)
        self.letter_order = (tuple(range(n_cls)) if letter_order is None else
                             tuple(letter_order))
        if sorted(self.letter_order) != list(range(n_cls)):
            raise ValueError(f"letter_order {self.letter_order} is not a "
                             f"permutation of the {n_cls} classes")
        self.reference = _load_reference()
        self.reference_strip = (config.VLM_REFERENCE_STRIP
                                if self.reference is not None else None)
        self.prompt = build_prompt(self.letter_order,
                                   self.reference_strip)
        self.budget = (config.VLM_TOKEN_BUDGET
                       if token_budget is None else token_budget)
        # padding_side="left" is load-bearing: with left padding the last
        # position of every row in a batch is the real end of that prompt, so
        # logits[:, -1] is the answer slot for all of them at once.
        self.processor = AutoProcessor.from_pretrained(self.model_id,
                                                       padding_side="left",
                                                       token=_hf_token())
        self._set_budget(self.budget)
        self.class_token_ids = self._resolve_class_tokens()
        self.n_image_tokens = self.measure_budget()

    def _set_budget(self, budget):
        """Pin the image processor's soft-token budget, by whichever name it
        uses. Records the names that took in `self.budget_attrs_set`; an
        empty list is not fatal on its own, because measure_budget is what
        actually decides whether the setting had an effect."""
        proc = getattr(self.processor, "image_processor", self.processor)
        self.budget_attrs_set = []
        for name in _BUDGET_ATTRS:
            if hasattr(proc, name):
                setattr(proc, name, budget)
                self.budget_attrs_set.append(name)

    def _resolve_class_tokens(self):
        """Token id for each class letter, asserting one token apiece.

        Both the bare letter and a leading-space variant are tried, since
        whether the answer slot starts a new word is tokenizer-dependent; the
        single-token spelling wins. Raises if a letter has none, because a
        multi-token label cannot be read out of one logit vector."""
        tok = self.processor.tokenizer
        ids = []
        for letter in config.VLM_CLASS_TOKENS:
            found = None
            for spelling in (letter, " " + letter):
                enc = tok.encode(spelling, add_special_tokens=False)
                if len(enc) == 1:
                    found = enc[0]
                    break
            if found is None:
                raise RuntimeError(
                    f"class label {letter!r} is not a single token for "
                    f"{self.model_id}; pick labels that are, or the "
                    "single-forward-pass readout cannot work")
            ids.append(found)
        if len(set(ids)) != len(ids):
            raise RuntimeError(f"class labels {config.VLM_CLASS_TOKENS} do "
                               "not map to distinct token ids")
        return ids

    def image_token_id(self):
        """The processor's image placeholder token id, or None.

        Tried processor-side first so the budget can be measured without a
        model; a caller holding one can pass its config as a fallback."""
        tid = getattr(self.processor, "image_token_id", None)
        if tid is not None:
            return int(tid)
        token = getattr(self.processor, "image_token", None)
        if token is not None:
            enc = self.processor.tokenizer.convert_tokens_to_ids(token)
            if enc is not None and enc >= 0:
                return int(enc)
        return None

    def measure_budget(self, model_config=None):
        """Count the image tokens the processor actually emitted for one
        probe image -- the honest answer to what the budget setting did.
        Returns -1 when the placeholder token cannot be identified."""
        probe = np.zeros((1, config.VLM_IMG_H, config.VLM_IMG_W, 3),
                         dtype=np.uint8)
        inputs = self.build_inputs(probe)
        tid = self.image_token_id()
        if tid is None and model_config is not None:
            tid = getattr(model_config, "image_token_id", None)
        if tid is None:
            return -1
        return int((inputs["input_ids"] == tid).sum())

    def report(self):
        """One-line summary of what the processor is actually doing, printed
        by every entry point so a surprise shows up before the run, not
        after."""
        if self.n_image_tokens < 0:
            print("  visual tokens per image: COULD NOT MEASURE (no image "
                  "token id on the processor) -- whether the requested "
                  f"budget of {self.budget} took is unverified.")
            return
        n_img = 2 if self.reference is not None else 1
        per = self.n_image_tokens / n_img
        print(f"  visual tokens: {self.n_image_tokens} over {n_img} image(s), "
              f"~{per:.0f} each (max_soft_tokens {self.budget}, set via "
              f"{config.VLM_BUDGET_KWARG or self.budget_attrs_set or 'nothing'})")
        # max_soft_tokens is a CEILING: the processor picks the largest patch
        # grid fitting under it, so landing below the request is correct and
        # only landing ABOVE it means the setting did not take.
        if per > self.budget + 1:
            print("  WARNING: more tokens per image than max_soft_tokens "
                  "allows -- the setting did not take, and prefill cost "
                  "follows the measured number.")

    def build_inputs(self, images):
        """Processor inputs for a batch of (n, H, W, 3) uint8 images.

        The arrays become PIL images in memory and go straight into the chat
        template -- there is no file and no base64 anywhere on this path.

        When a reference figure is configured it is the FIRST image in every
        conversation and the event is the last, so `Classify the single event
        in the last image` is unambiguous. That reference is a second image
        on every request, which is not free: it is paid on each forward pass
        because a plain batched forward has no prefix cache. Since it is
        byte-identical and leads the prompt it is a genuine shared prefix, so
        if it earns its keep the next optimization is caching its KV once and
        expanding across the batch.
        """
        from PIL import Image

        messages = []
        for im in images:
            content = []
            # Order is reference, then TEXT, then the event -- deliberately
            # not the model card's images-first layout. Everything before the
            # event image is byte-identical on every request, so this makes
            # the reference plus the whole instruction one shared prefix that
            # can be cached (see EventClassifier._build_prefix_cache); with
            # the text last, only the reference would be shared and the cache
            # would save a fraction as much.
            if self.reference is not None:
                content.append({"type": "image", "image": self.reference})
            content.append({"type": "text", "text": self.prompt})
            content.append({"type": "image", "image": Image.fromarray(im)})
            messages.append([{"role": "user", "content": content}])
        proc_kwargs = {"return_tensors": "pt", "padding": True}
        if config.VLM_BUDGET_KWARG:
            proc_kwargs[config.VLM_BUDGET_KWARG] = self.budget
        return self.processor.apply_chat_template(
            messages,
            add_generation_prompt=True,
            tokenize=True,
            return_dict=True,
            processor_kwargs=proc_kwargs)


class EventClassifier:
    """A VLM held open across batches, answering with four class posteriors.

    Construction loads the model behind a PromptBuilder, resolves the kwarg
    that keeps the LM head off every position but the answer slot, and
    reports what the processor actually produced.
    """

    def __init__(self, model_id=None, token_budget=None, letter_order=None,
                 builder=None):
        import torch

        self._torch = torch
        self.builder = builder or PromptBuilder(model_id, token_budget,
                                                letter_order)
        self.model_id = self.builder.model_id
        self.letter_order = self.builder.letter_order
        self.class_token_ids = self.builder.class_token_ids
        print(f"vlm: loading {self.model_id} "
              f"({config.VLM_DTYPE}, {config.VLM_DEVICE})")
        self.model = self._load_model()
        self.model.eval()
        self._logits_kwarg = self._resolve_logits_kwarg()
        if self.builder.n_image_tokens < 0:
            self.builder.n_image_tokens = self.builder.measure_budget(
                self.model.config)
        self.builder.report()
        self._prefix = None
        self._prefix_batch = -1
        self._prefix_unverified = None
        self._timing = {"render": 0.0, "processor": 0.0, "forward": 0.0}
        self._timed = 0
        if config.VLM_PREFIX_CACHE:
            self._setup_prefix_cache(config.VLM_BATCH)

    def _setup_prefix_cache(self, batch):
        """Build and verify the shared-prefix cache for one batch width.

        Runs before any classification, so an OOM here must NOT kill the run:
        the cache is an optimization, and without it classify() still works
        and still backs off. Failing loudly but continuing is the point.
        """
        torch = self._torch
        try:
            state = self._build_prefix_cache(batch)
            if state is None:
                return
            # Verified on the first real batch, not here: synthetic probes
            # cannot tell us whether real tags move.
            self._prefix_unverified = state
        except torch.OutOfMemoryError:
            torch.cuda.empty_cache()
            print(f"  prefix cache: OOM building it at batch {batch}; "
                  "continuing uncached. classify() will retry at whatever "
                  "batch it backs off to.")
            self._prefix = None
            self._prefix_batch = -1
            self._prefix_unverified = None
            return
        except (RuntimeError, ValueError, KeyError) as exc:
            # A broken cache must not cost a run. The reason is printed in
            # full rather than swallowed, because this path failing means
            # the prefix split is wrong and wants fixing, not ignoring.
            torch.cuda.empty_cache()
            print(f"  prefix cache: DISABLED, {type(exc).__name__}: {exc}")
            self._prefix = None
            self._prefix_batch = -1
            self._prefix_unverified = None
            return
        self._prefix_batch = batch
        print(f"  prefix cache: built at batch {batch}, "
              f"{state['start']} tokens shared per event; verifying on the "
              "first real batch")

    def _load_model(self):
        """Load the checkpoint, quantized if config.VLM_LOAD_IN asks.

        The auto class differs across the ladder: the E-series cards name
        AutoModelForImageTextToText, while 12B is the Unified (encoder-free)
        variant and names AutoModelForMultimodalLM. Both are tried rather
        than hard-coding one, so moving a rung does not need a code change.
        """
        import torch
        import transformers

        kw = {
            "device_map": config.VLM_DEVICE,
            "attn_implementation": config.VLM_ATTN,
            "token": _hf_token(),
        }
        want = getattr(config, "VLM_LOAD_IN", None)
        if want:
            from transformers import BitsAndBytesConfig
            if want == "8bit":
                kw["quantization_config"] = BitsAndBytesConfig(
                    load_in_8bit=True)
            elif want == "4bit":
                kw["quantization_config"] = BitsAndBytesConfig(
                    load_in_4bit=True,
                    bnb_4bit_quant_type="nf4",
                    bnb_4bit_compute_dtype=getattr(torch, config.VLM_DTYPE))
            else:
                raise ValueError(f"VLM_LOAD_IN={want!r}; want '8bit', "
                                 "'4bit' or None")
            # device_map must let accelerate place the quantized shards.
            kw["device_map"] = "auto"
        else:
            kw["dtype"] = getattr(torch, config.VLM_DTYPE)

        last = None
        for name in ("AutoModelForImageTextToText",
                     "AutoModelForMultimodalLM"):
            cls = getattr(transformers, name, None)
            if cls is None:
                continue
            try:
                model = cls.from_pretrained(self.model_id, **kw)
                print(f"  loaded via {name}"
                      + (f", {want}" if want else ", bf16"))
                return model
            except (ValueError, KeyError, TypeError) as exc:
                last = exc
        raise RuntimeError(
            f"no auto class could load {self.model_id}: {last}")

    def _report_timing(self):
        """Where the wall clock actually goes, as a share of the three
        stages. `render` is the numpy rasterizer, `processor` is the chat
        template plus Gemma's image preprocessing, `forward` is the GPU."""
        total = sum(self._timing.values())
        if total <= 0.0 or self._timed == 0:
            return
        parts = "  ".join(f"{k} {100.0 * v / total:4.1f}%"
                          for k, v in self._timing.items())
        print(f"  timing: {parts}  |  {self._timed / total:5.1f} ev/s")

    def _resolve_logits_kwarg(self):
        """Name of the forward kwarg that limits the LM head to the last
        positions, or None if this model takes neither spelling.

        transformers renamed `num_logits_to_keep` to `logits_to_keep`; both
        are in the wild. Returning None is survivable but expensive -- see
        class_probs -- so it is reported rather than passed over.
        """
        import inspect
        try:
            params = inspect.signature(self.model.forward).parameters
        except (TypeError, ValueError):
            return None
        for name in ("logits_to_keep", "num_logits_to_keep"):
            if name in params:
                return name
        print("  WARNING: this model takes neither logits_to_keep nor "
              "num_logits_to_keep, so the LM head runs over every position. "
              "Drop config.VLM_BATCH if the forward pass runs out of memory.")
        return None

    # --- shared-prefix KV cache ------------------------------------------
    #
    # build_inputs orders every conversation as [reference][text][event], so
    # every token before the event's own image block is identical on each
    # request. Encoding that once and reusing its keys and values leaves only
    # the event's ~62 image tokens to compute per event.
    #
    # The cache is built at a fixed batch size and cropped back to the prefix
    # length after each use, rather than having its tensors repeated by hand:
    # a Cache is an object with a layout that changes between transformers
    # releases, and crop() is the supported way to rewind one. A batch that
    # does not match the cached width falls back to the plain path.

    def _image_token_runs(self, ids):
        """Contiguous runs of image-placeholder tokens in one row, as
        (start, stop) pairs. There are two when a reference figure is in
        play -- the reference, then the event -- and one when it is not."""
        tid = self.builder.image_token_id()
        if tid is None:
            return []
        runs = []
        for i in (ids == tid).nonzero().flatten().tolist():
            if runs and i == runs[-1][1]:
                runs[-1][1] = i + 1
            else:
                runs.append([i, i + 1])
        return [tuple(r) for r in runs]

    def _split_inputs(self, inputs, start):
        """Per-TOKEN tensors sliced from `start`; everything else passed
        through. Image-indexed tensors are dropped here and re-added by the
        caller via _image_slice, because they are indexed by image rather
        than by token and slicing them by token silently mismatches the two.
        """
        out = {}
        n_tok = inputs["input_ids"].shape[1]
        for k, v in inputs.items():
            if k in _IMAGE_INDEXED:
                continue
            if hasattr(v, "ndim") and v.ndim >= 2 and v.shape[1] == n_tok:
                out[k] = v[:, start:]
            else:
                out[k] = v
        return out

    def _image_slice(self, inputs, batch, first=False):
        """The event images' image-indexed tensors, or the reference's with
        first=True.

        Images are stacked in conversation order, so with a reference each
        sample contributes [reference, event] and the event is every second
        entry. Both pixel_values and image_position_ids are indexed this way:
        Gemma 4's patch embedder adds a position embedding per image, so
        handing it 32 images' pixels alongside 64 images' positions fails on
        a shape mismatch inside the vision tower.
        """
        per = 2 if self.builder.reference is not None else 1
        out = {}
        for k in _IMAGE_INDEXED:
            v = inputs.get(k)
            if v is None:
                continue
            if v.shape[0] != batch * per:
                raise RuntimeError(
                    f"{k} has {v.shape[0]} entries for {batch} events at "
                    f"{per} images each; cannot tell the reference from the "
                    "event")
            out[k] = v[0::per] if first else v[per - 1::per]
        return out

    def _build_prefix_cache(self, batch):
        """Encode the shared prefix once at this batch width, or None."""
        torch = self._torch
        probe = np.zeros((batch, config.VLM_IMG_H, config.VLM_IMG_W, 3),
                         dtype=np.uint8)
        inputs = self.builder.build_inputs(probe).to(self.model.device)
        runs = self._image_token_runs(inputs["input_ids"][0])
        if not runs:
            print("  prefix cache: no image tokens found, staying uncached")
            return None
        start = runs[-1][0]  # the event's image block always comes last
        if start <= 0:
            print("  prefix cache: nothing precedes the event, nothing to "
                  "cache")
            return None
        head = {}
        n_tok = inputs["input_ids"].shape[1]
        for k, v in inputs.items():
            if k in _IMAGE_INDEXED:
                continue
            if hasattr(v, "ndim") and v.ndim >= 2 and v.shape[1] == n_tok:
                head[k] = v[:, :start]
            else:
                head[k] = v
        head.update(self._image_slice(inputs, batch, first=True))
        # Build the cache from the model config rather than taking whatever
        # the forward defaults to. Gemma 4 interleaves local sliding-window
        # attention with global attention, and a cache that does not know
        # which layers are which caches the wrong span for the local ones --
        # which is the shape of a large logit difference on an otherwise
        # correct split.
        cache = None
        try:
            from transformers import DynamicCache
            cache = DynamicCache(config=self.model.config)
        except (ImportError, TypeError) as exc:
            print(f"  prefix cache: no config-aware DynamicCache ({exc}); "
                  "falling back to the forward's default cache")
        with torch.inference_mode():
            out = self.model(**head, use_cache=True,
                             **({"past_key_values": cache}
                                if cache is not None else {}))
        cache = out.past_key_values
        if not hasattr(cache, "crop"):
            print("  prefix cache: this Cache has no crop(), staying "
                  "uncached")
            return None
        return {"cache": cache, "start": start, "batch": batch}

    def _verify_prefix_cache(self, state, images):
        """Accept the cache only if it changes no DECISION, on real events.

        Two earlier criteria were wrong and worth recording. Comparing raw
        logits rejects a usable cache: logits run to order +-20 and bfloat16
        accumulating over ~240 tokens in one pass differs from ~62 with a
        cache in the last digits (~0.8), which the softmax washes out.
        Comparing all four posteriors on synthetic all-black / all-white
        probes is barely better -- it scores classes nothing depends on, at
        pixel values no event ever has.

        What the yield rests on is the (a,n) tag at
        config.VLM_AN_THRESHOLD, so that is the test: on a real batch, no
        event may change its tag, and |dp(an)| must stay under tolerance.
        The images must also produce differing posteriors, or a cache that
        stopped routing the event image through would pass trivially.
        """
        plain = self._posteriors(self._forward_logits(images, state=None))
        cached = self._posteriors(self._forward_logits(images, state=state))
        d_an = float(np.abs(an_probability(plain)
                            - an_probability(cached)).max())
        flips = int((is_an(plain) != is_an(cached)).sum())
        pan = an_probability(plain)
        spread = float(np.abs(pan - pan.mean()).max())
        # How much of this batch could plausibly have flipped: events sitting
        # within the measured shift of the threshold. Zero flips out of a
        # batch where nothing is near the boundary is not evidence, and
        # saying so is more useful than a tolerance constant.
        near = int((np.abs(pan - config.VLM_AN_THRESHOLD) <= d_an).sum())
        print(f"  prefix cache: on {images.shape[0]} real events, "
              f"tag flips {flips}, max |dp(an)| {d_an:.4f}, "
              f"{near} event(s) within that of the p>={config.VLM_AN_THRESHOLD}"
              " threshold")
        if spread < 1e-4:
            print("  prefix cache: DISABLED -- every event scores the same "
                  "p(an), so this batch cannot test anything.")
            return False
        if flips:
            print("  prefix cache: DISABLED -- it changes tags, so it is a "
                  "bug rather than an optimization.")
            return False
        if near == 0:
            print(f"  prefix cache: live, but NOTE no event in this batch "
                  f"was within {d_an:.3f} of the threshold, so zero flips is "
                  "weak evidence. The shift is real; it only costs anything "
                  "for events near the cut.")
        else:
            print("  prefix cache: verified, live")
        if d_an > config.VLM_PREFIX_CACHE_TOL:
            print(f"  prefix cache: p(an) moves up to {d_an:.3f}, above the "
                  f"{config.VLM_PREFIX_CACHE_TOL} noted in config. That is a "
                  "systematic on any threshold near the bulk of the p(an) "
                  "distribution -- compare it against the seed spread before "
                  "quoting either.")
        return True

    def _forward_logits(self, images, state):
        """Logits at the answer slot, with or without the prefix cache.

        The timing wraps BOTH paths. An earlier version incremented the
        forward accumulator and the event counter only in the uncached
        branch, so switching the cache on froze the denominator and made a
        working cache report a collapsing event rate -- the instrumentation
        decaying, not the run.
        """
        torch = self._torch
        timed = bool(config.VLM_TIMING_EVERY)
        t0 = time.perf_counter() if timed else 0.0
        inputs = self.builder.build_inputs(images).to(self.model.device)
        if timed:
            torch.cuda.synchronize()
            self._timing["processor"] += time.perf_counter() - t0
            t0 = time.perf_counter()
        kw = {}
        if self._logits_kwarg:
            kw[self._logits_kwarg] = 1
        try:
            if state is None:
                with torch.inference_mode():
                    return self.model(**inputs,
                                      **kw).logits[:, -1, :].float()
            cut, cache = state["start"], state["cache"]
            tail = self._split_inputs(inputs, cut)
            tail.update(self._image_slice(inputs, images.shape[0]))
            # attention_mask stays FULL length: the cached prefix is still
            # being attended to even though its tokens are not in this
            # forward.
            tail["attention_mask"] = inputs["attention_mask"]
            tail["past_key_values"] = cache
            try:
                with torch.inference_mode():
                    out = self.model(**tail, use_cache=True, **kw)
                return out.logits[:, -1, :].float()
            finally:
                # Rewind so the next batch reuses the prefix. A sliding-window
                # layer refuses to be cropped once it has seen more tokens
                # than its window (512 for Gemma 4), so this works only while
                # prefix + suffix stays under it -- at max_soft_tokens 70 that
                # is 389 + 62, at 140 it is 455 + 120 and it does not. The
                # forward already completed, so this batch's result stands;
                # the cache is simply dropped for every batch after it.
                try:
                    cache.crop(cut)
                except Exception as exc:
                    print(f"  prefix cache: DISABLED after use, "
                          f"{type(exc).__name__}: {exc}")
                    self._prefix = None
                    self._prefix_batch = -1
                    self._prefix_unverified = None
        finally:
            if timed:
                torch.cuda.synchronize()
                self._timing["forward"] += time.perf_counter() - t0
                self._timed += images.shape[0]

    def class_probs(self, images):
        """Class posteriors (n, len(VLM_CLASSES)) float32 for a batch.

        One forward pass; the logits at the answer slot are gathered at the
        four class token ids and softmaxed over those four alone, so the
        probabilities are conditional on the answer being one of the classes.
        Column order follows config.VLM_CLASSES.

        `logits_to_keep=1` is not an optimization, it is what makes the batch
        fit: Gemma's vocabulary is ~262k, so running the LM head over every
        position would allocate batch x seq x 262144 -- gigabytes of logits
        per forward, on a card that has to hold the model too. Only the
        answer slot is ever read.
        """
        if (self._prefix_unverified is not None
                and images.shape[0] == self._prefix_batch):
            pending, self._prefix_unverified = self._prefix_unverified, None
            try:
                ok = self._verify_prefix_cache(pending, images)
            except Exception as exc:
                print(f"  prefix cache: DISABLED, {type(exc).__name__}: "
                      f"{exc}")
                ok = False
            if ok and self._prefix_batch > 0:
                self._prefix = pending
            else:
                self._prefix = None
                self._prefix_batch = -1
        state = (self._prefix
                 if images.shape[0] == self._prefix_batch else None)
        return self._posteriors(self._forward_logits(images, state))

    def _posteriors(self, logits):
        """Logits at the answer slot -> class posteriors in VLM_CLASSES order.

        class_token_ids are in LETTER order; softmax over the letters, then
        scatter back to canonical order so every caller sees the same columns
        whatever permutation ran."""
        torch = self._torch
        picked = logits[:, self.class_token_ids]
        probs = torch.softmax(picked, dim=-1).cpu().numpy()
        out = np.empty_like(probs)
        out[:, list(self.letter_order)] = probs
        return out

    def classify(self, X, beam_ref=None, X2=None, batch=None, dy=0.0,
                 dhalf=0.0):
        """Render and classify a full set of traces, batch by batch.

        Returns (probs (n, n_classes) float32, labels (n,) int64), where the
        label is the argmax column. Only the (a,n) column decides anything
        downstream -- see an_probability / is_an -- but the full posterior is
        returned so a threshold can be swept without re-running.

        On CUDA OOM the batch is halved and the chunk retried, permanently,
        rather than the run dying. Gemma 4's vision tower builds a one_hot
        over pre-pool patch positions, which at max_soft_tokens=280 is a
        48x48 grid per image and doubles again when a reference figure rides
        along -- so the workable batch is much smaller than the token count
        alone suggests, and is easier to discover than to predict.
        """
        torch = self._torch
        n = X.shape[0]
        size = batch or config.VLM_BATCH
        if config.VLM_PREFIX_CACHE and size != self._prefix_batch:
            self._prefix = None
            self._prefix_batch = -1
            self._setup_prefix_cache(size)
        out = np.empty((n, len(config.VLM_CLASSES)), dtype=np.float32)
        start = 0
        batches = 0
        while start < n:
            stop = min(start + size, n)
            try:
                chunk2 = None if X2 is None else X2[start:stop]
                t0 = time.perf_counter()
                images = render_traces(X[start:stop], beam_ref, chunk2, dy,
                                       dhalf)
                if config.VLM_TIMING_EVERY:
                    self._timing["render"] += time.perf_counter() - t0
                out[start:stop] = self.class_probs(images)
            except torch.OutOfMemoryError:
                torch.cuda.empty_cache()
                if size == 1:
                    raise RuntimeError(
                        "CUDA OOM at batch size 1. Lower "
                        "config.VLM_TOKEN_BUDGET (the tiers are 70, 140, "
                        "280, 560, 1120) -- memory in the vision tower grows "
                        "with the square of the patch grid, so one tier down "
                        "is roughly a quarter the cost.") from None
                size = max(1, size // 2)
                print(f"  CUDA OOM, batch -> {size}          ")
                if self._prefix is not None:
                    # The cache is built at one width; rebuild it at the new
                    # one or every batch silently falls back to plain.
                    self._prefix = None
                    self._prefix_batch = -1
                    self._setup_prefix_cache(size)
                continue
            start = stop
            batches += 1
            if config.VLM_TIMING_EVERY and \
                    batches % config.VLM_TIMING_EVERY == 0:
                print(f"  classified {stop}/{n}")
                self._report_timing()
            else:
                print(f"  classified {stop}/{n}", end="\r", flush=True)
        print(f"  classified {n}/{n}      ")
        return out, out.argmax(axis=1).astype(np.int64)
