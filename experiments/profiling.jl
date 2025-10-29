using Profile, Profile.Allocs, PProf
# open the allocation UI from current session data
PProf.Allocs.pprof(web=true, webhost="127.0.0.1", webport=0)