let
  ethan-desktop = "ssh-ed25519 AAAAC3NzaC1lZDI1NTE5AAAAIDlbs+h9OqZMIAC6b3i4tUcXC4PidfBFEQNdwrLS8g9G";
  allKeys = [
    ethan-desktop
  ];
in
{
  "settings.json.age".publicKeys = allKeys;
}
