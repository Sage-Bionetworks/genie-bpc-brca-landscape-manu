# This is not actually a markdown file.  I just did that so it didn't get ignored.

dag {
  bb="-2.929,-4.493,2.771,3.689"
  Age [selected,pos="-0.940,-0.313"]
  Gene_feature [exposure,pos="-0.098,1.042"]
  OS [outcome,pos="0.922,1.052"]
  Sample_type [pos="-1.226,1.024"]
  Unmeas_confounder [pos="-0.903,2.501"]
  Age -> Gene_feature
  Age -> OS
  Gene_feature -> OS
  Sample_type -> Gene_feature
  Unmeas_confounder -> Gene_feature
  Unmeas_confounder -> OS
}
