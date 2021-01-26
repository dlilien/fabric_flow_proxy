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
shift_result.sh rstf/d_ftw_rc2.0e-3_oop0.0_bm1.0e-2_0013.vtu rstf/d_ftw_rc2.0e-3_oop0.0_bm1.0e-2_r_????.vtu

# shift_result.sh rstf/d_ftw_rc2.0e-3_oop0.0_bm1.0e-2_0015.vtu rstf/d_ftw_rc2.0e-3_oop0.0_bm1.0e-2_rr_????.vtu
touch rstf/d_ftw_rc2.0e-3_oop0.0_bm1.0e-2.result

shift_result.sh rstf/d_ftw_rc2.0e-2_oop0.0_bm0.0_0004.vtu rstf/d_ftw_rc2.0e-2_oop0.0_bm0.0_r_????.vtu
touch rstf/d_ftw_rc2.0e-2_oop0.0_bm0.0.result


# shift_result.sh rstf/d_ftw_rc2.0e-4_oop0.0_bm0.0_2315.vtu rstf/d_ftw_rc2.0e-4_oop0.0_bm0.0_r_????.vtu

shift_result.sh stream_ftw/d_ftw_rc2.0e-3_oop1.0e-3_bm0.0_0003.vtu  stream_ftw/d_ftw_rc2.0e-3_oop1.0e-3_bm0.0_r_????.vtu
touch stream_ftw/d_ftw_rc2.0e-3_oop1.0e-3_bm0.0.result
