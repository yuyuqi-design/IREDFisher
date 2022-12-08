change HIE to HIS
change ILE HG2,HG3 to HG1,HG2  sed -i 's/2HG2 ILE/1HG2 ILE/;s/3HG2 ILE/2HG2 ILE' file.pdb




python2 /home/g02808yy/software/rosetta2019/main/source/scripts/python/public/molfile_to_params.py AMP.mol2 -p AMP -n AMP --keep-names --clobber
python make_flags.py AMB_CarNi.pdb 297A,507A,613A,637B,638C 609A,611A,392A,389A,338A,299A,418A,277A,423A,426A
sh make_csts.sh AMB_CarNi.pdb > bbCA.cst
mpirun.mpich -np 20 rosetta_scripts.mpi.linuxgccrelease @refine.flags -nstruct 20
for i in `seq 0 1125 22500`;do head -$i jobs.sh | tail -1125 > $i'_job.sh';done
grep AMB_relaxed score_0*|sort -n -k2|head



residues library /home/g02808yy/Downloads/rosetta_src_2020.08.61146_bundle/main/database/chemical/residue_type_sets/fa_standard/residue_types/
