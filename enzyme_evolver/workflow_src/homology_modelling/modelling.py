from enzyme_evolver.workflow_src.homology_modelling import auto_modeller

class model_class:
    def __init__(self, seqcode='', db='enzyme_evolver/database/pdb_95.pir', auto_single='auto_single', auto_multiple='auto_multiple', template_single='template_single', homodimer='multi_chain', start_residue='start_residue',end_residue='end_residue'):
        self.seqcode = seqcode
        self.db = db
        self.auto_single = auto_single
        self.auto_multiple = auto_multiple
        self.template_single = template_single
        self.homodimer =homodimer
        self.start_residue=start_residue
        self.end_residue=end_residue

    def auto_single_modelling(self):
        auto_modeller.template_search(self.seqcode,self.db)
        templates, template1, template2, template3, template1_chain, template1_PdbCode, template2_chain, template2_PdbCode, template3_chain, template3_PdbCode = auto_modeller.check_template(self.seqcode)
        auto_modeller.align2d(self.seqcode,template1,template1_PdbCode,template1_chain)
        auto_modeller.build_model(self.seqcode,template1)

    def auto_multiple_modelling(self):
        auto_modeller.template_search(self.seqcode,self.db)
        templates, template1, template2, template3, template1_chain, template1_PdbCode, template2_chain, template2_PdbCode, template3_chain, template3_PdbCode = auto_modeller.check_template(self.seqcode)
        auto_modeller.salign(self.seqcode, template1_PdbCode,template1_chain,template2_PdbCode, template2_chain,template3_PdbCode, template3_chain)
        auto_modeller.align2d_mult(self.seqcode)
        auto_modeller.build_multi_models(self.seqcode,template1,template2,template3)

    def template_single_modelling(self):
        template1 = self.template_single.split('.')[0] + 'A'
        template1_PdbCode = self.template_single.split('.')[0]
        template1_chain = 'A'
        auto_modeller.align2d(self.seqcode,template1,template1_PdbCode,template1_chain)
        auto_modeller.build_model(self.seqcode,template1)
    def homodimer_modelling(self):
        auto_modeller.single2dimer_ali(self.seqcode, self.start_residue,self.end_residue,'2')
        auto_modeller.build_chains(self.seqcode)