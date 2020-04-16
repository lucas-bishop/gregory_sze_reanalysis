mothur "#get.groups(shared=final.shared, groups=z_mock_ext-Zymo_mockpcr_1-Zymo_mockpcr_2-Zymo_mockpcr_3)"

# Renaming output file
mv final.0.03.pick.shared mock.final.shared



# Control shared file
echo PROGRESS: Creating control shared file.

# Removing any non-control groups from shared file
mothur "#get.groups(shared=final.shared, groups=PCRwater_neg_1-PCRwater_neg_2-PCRwater_neg_3)"

# Renaming output file
mv final.0.03.pick.shared control.final.shared


