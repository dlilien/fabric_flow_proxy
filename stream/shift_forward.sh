#! /bin/bash
#
# shift_bm.sh
# Copyright (C) 2020 dal22 <dal22@fenrir.ess.washington.edu>
#
# Distributed under terms of the MIT license.
#

# Shift results with centimeters of basal melt since I had to restart due to weird sensitivity after 100-200 years.

# We intentionally ignore 2 so we don't have problems with the initial steps we are dumping
# This should just cause a linking error, and then we continue on our merry way
# shift_result.sh rstf/s_ftw_rc2.0e-3_oop1.0e-1_bm0.0_0025.vtu rstf/s_ftw_rc2.0e-3_oop1.0e-1_bm0.0_r_????.vtu
# shift_result.sh rstf/s_ftw_rc2.0e-3_oop1.0e-1_bm0.0_0102.vtu rstf/s_ftw_rc2.0e-3_oop1.0e-1_bm0.0_rr_????.vtu
# touch rstf/s_ftw_rc2.0e-3_oop1.0e-1_bm0.0.result


# shift_result.sh rstf/s_nuc_rc2.0e-3_oop0.0_bm0.0_0947.vtu rstf/s_nuc_rc2.0e-3_oop0.0_bm0.0_r_????.vtu
# touch rstf/s_nuc_rc2.0e-3_oop0.0_bm0.0.result

# shift_result.sh rstf/s_ftw_rc2.0_oop0.0_bm0.0_0953.vtu rstf/s_ftw_rc2.0_oop0.0_bm0.0_r_????.vtu
# touch rstf/s_ftw_rc2.0_oop0.0_bm0.0.result

shift_result.sh rstf/s_ftw_rc2.0e-3_oop1.0e-1_bm0.0_0237.vtu rstf/s_ftw_rc2.0e-3_oop1.0e-1_bm0.0_rrr_????.vtu
touch rstf/s_ftw_rc2.0e-3_oop1.0e-1_bm0.0.result

exit

# shift_result.sh rstf/s_ftw_rc2.0e-3_oop0.0_bm1.0e-3_0907.vtu rstf/s_ftw_rc2.0e-3_oop0.0_bm1.0e-3_r_????.vtu
# touch rstf/s_ftw_rc2.0e-3_oop0.0_bm1.0e-3.result

# shift_result.sh rstf/s_ftw_rc2.0e-3_oop0.0_bm1.0e-2_0887.vtu rstf/s_ftw_rc2.0e-3_oop0.0_bm1.0e-2_r_????.vtu
# touch rstf/s_ftw_rc2.0e-3_oop0.0_bm1.0e-2.result

shift_result.sh rstf/s_ftw_rc2.0e-3_oop0.0_bm2.0e-2_0882.vtu rstf/s_ftw_rc2.0e-3_oop0.0_bm2.0e-2_r_????.vtu
touch rstf/s_ftw_rc2.0e-3_oop0.0_bm2.0e-2.result

# shift_result.sh rstf/s_ftw_rc2.0e-3_oop1.0e-2_bm0.0_0373.vtu rstf/s_ftw_rc2.0e-3_oop1.0e-2_bm0.0_r_????.vtu
# touch rstf/s_ftw_rc2.0e-3_oop1.0e-2_bm0.0.result

# shift_result.sh rstf/s_ftw_rc2.0e-3_oop1.0e-3_bm0.0_0314.vtu rstf/s_ftw_rc2.0e-3_oop1.0e-3_bm0.0_r_????.vtu
# touch rstf/s_ftw_rc2.0e-3_oop1.0e-3_bm0.0.result

# shift_result.sh rstf/s_nuclc_rc2.0e-3_oop0.0_bm0.0_0650.vtu rstf/s_nuclc_rc2.0e-3_oop0.0_bm0.0_r_????.vtu
# touch rstf/s_nuclc_rc2.0e-3_oop0.0_bm0.0.result

# shift_result.sh rstf/s_ftwlc_rc2.0e-3_oop0.0_bm0.0_0952.vtu rstf/s_ftwlc_rc2.0e-3_oop0.0_bm0.0_r_????.vtu
# touch rstf/s_ftwlc_rc2.0e-3_oop0.0_bm0.0.result

# shift_result.sh rstf/s_ftw_rc2.0e-3_oop0.0_bm0.0_0923.vtu rstf/s_ftw_rc2.0e-3_oop0.0_bm0.0_r_????.vtu
# touch rstf/s_ftw_rc2.0e-3_oop0.0_bm0.0.result

# shift_result.sh rstf/s_ftw_rc2.0e-4_oop0.0_bm0.0_0914.vtu rstf/s_ftw_rc2.0e-4_oop0.0_bm0.0_r_????.vtu
# touch rstf/s_ftw_rc2.0e-4_oop0.0_bm0.0.result

# shift_result.sh rstf/s_ftw_rc2.0e-2_oop0.0_bm0.0_0095.vtu rstf/s_ftw_rc2.0e-2_oop0.0_bm0.0_r_????.vtu
# touch rstf/s_ftw_rc2.0e-2_oop0.0_bm0.0.result

# shift_result.sh rstf/s_ftw_rc2.0e-1_oop0.0_bm0.0_0876.vtu rstf/s_ftw_rc2.0e-1_oop0.0_bm0.0_r_????.vtu
# touch rstf/s_ftw_rc2.0e-1_oop0.0_bm0.0.result
