"""Writing TagEfficiencyRecord stores that the C++ cross section can read.

TagEfficiency.hpp states the contract plainly: the efficiency side hands the
cross section one record per channel and reaction strip -- the count, the
fraction of true reactions that end up in it, and the fraction counted one
strip late -- and "change the method, keep the record, and the cross section
keeps working". So a Python classifier does not need a new path into
CrossSection; it needs to write this file.

The on-disk format is whatever TagEfficiencyStore::Write produces, mirrored
here exactly: a flat ROOT file of TParameter<double> objects keyed
`<channel>_<field>_r<reac>` over the five fields, plus a TNamed "method"
stamped so a later reader knows what it is holding. Records for the channel
being written are deleted first, so a strip that lost its tag does not keep a
stale one.

Which file the C++ reads is CROSS_SECTION_CONFIG.TAG_EFFICIENCY_FILE (a
basename under root_files/, defaulting to tag_efficiency.root). Writing the
VLM records to their own basename therefore leaves the bootstrap store
untouched and makes the two methods a one-line config switch apart -- which
is also how a per-seed ensemble is run.
"""

import config

_FIELDS = ("n_counted", "eff", "eff_err", "migrate", "tag_eff")
_K_OVERWRITE = 2  # TObject::kOverwrite


def store_path(basename):
    """Absolute path of a store, given its basename under root_files/."""
    return str(config.ROOT_FILES_DIR / basename)


def write(basename, channel, records, method):
    """Write one channel's records, mirroring TagEfficiencyStore::Write.

    `records` is an iterable of dicts with a `reac` key plus the five record
    fields; a missing field is written as 0.0, which is what an unset
    TagEfficiencyRecord field holds on the C++ side. Opens UPDATE, so
    channels already in the file survive.
    """
    import ROOT

    path = store_path(basename)
    f = ROOT.TFile(path, "UPDATE")
    if f.IsZombie():
        raise RuntimeError(f"cannot write tag-efficiency store {path}")
    f.cd()
    ROOT.TNamed("method", method).Write("method", _K_OVERWRITE)
    # Drop this channel's previous records before writing the new ones, so a
    # strip that is no longer tagged does not keep a stale record.
    prefix = f"{channel}_"
    for key in [k.GetName() for k in f.GetListOfKeys()]:
        if key.startswith(prefix):
            f.Delete(f"{key};*")
    n = 0
    for rec in records:
        reac = int(rec["reac"])
        for field in _FIELDS:
            key = f"{channel}_{field}_r{reac}"
            param = ROOT.TParameter(float)(key, float(rec.get(field, 0.0)))
            param.Write(key, _K_OVERWRITE)
        n += 1
    f.Close()
    print(f"  wrote {n} {channel} records to {path} (method: {method})")
    return path
