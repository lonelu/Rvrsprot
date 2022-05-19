To generate C2 4 helix bundles that bind a porphyrin.

1. run 'gss_sc_rots.py' to get the rotations. Result in directory '../std_rots/'.
2. run 'gss_zfix_mvdmH.py' to get 1 helix bind the metal with z angle fixed.
3. run 'gss_zfix_vdmH.py' to get the other helix form a 2nd shell with z angle fixed.
4. run 'gss_c2.py' to combine the two helixs and rotate to form 4 helix bundles.
5. run cccp fit in cccp matlab program. (a. generate_CA_coords.py; b. run_fit.m; c. _fcick_read.py.)
6. run 'extract_good_cccp.py' to copy the structures with low rmsd with cccp models.