#! /bin/bash
#
# sync_to_fenrir.sh
# Copyright (C) 2020 dal22 <dal22@loki.ess.washington.edu>
#
# Distributed under terms of the MIT license.
#


# rsync stream_ftw/s_*vtu stream_ftw/s_*result /colddata/loki-rd04/dal22/anisotropy/stream/stream_ftw/
rsync -L stream_ftw/*.tar.gz /colddata/loki-rd04/dal22/anisotropy/stream/stream_ftw/
rsync -L rstf/*.tar.gz /colddata/loki-rd04/dal22/anisotropy/stream/rstf/
