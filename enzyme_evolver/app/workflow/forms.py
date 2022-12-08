from flask_wtf import FlaskForm, Form
from flask_wtf.file import FileField, FileRequired, FileAllowed
from wtforms import StringField, BooleanField, SubmitField, IntegerField, DecimalField, SelectField, RadioField, TextAreaField, MultipleFileField,FormField, RadioField
from wtforms.validators import DataRequired, NumberRange, ValidationError,regexp, Optional




AA_LIST = ['X', 'V', 'G', 'F', 'E', 'N', 'P', 'Q', 'M', 'K', 'T', 'S', 'W', 'A', 'R', 'D', 'L', 'Y', 'H',
                    'I', 'C', '*',' ','\n','\r']

def is_protein_sequence(form, field):
    for aa in field.data.upper():
        if aa not in AA_LIST:
            raise ValidationError('Non-amino acid characters are not accepted')



class RequiredIf(DataRequired):

    """Validator which makes a field required if another field is set and has a truthy value.

    Sources:
        - http://wtforms.simplecodes.com/docs/1.0.1/validators.html
        - http://stackoverflow.com/questions/8463209/how-to-make-a-field-conditionally-optional-in-wtforms
        - https://gist.github.com/devxoul/7638142#file-wtf_required_if-py
    """
    field_flags = ('requiredif',)

    def __init__(self, message=None, *args, **kwargs):
        super(RequiredIf).__init__()
        self.message = message
        self.conditions = kwargs

    # field is requiring that name field in the form is data value in the form
    def __call__(self, form, field):
        for name, data in self.conditions.items():
            other_field = form[name]
            if other_field is None:
                raise Exception('no field named "%s" in form' % name)
            if other_field.data == data and not field.data:
                DataRequired.__call__(self, form, field)
            Optional()(form, field)

class ExampleForm(Form):
    example_check = BooleanField('Load example')

class NewWorkFlowForm(FlaskForm):
    name = StringField(validators=[DataRequired()])
    submit = SubmitField('Submit')

class HomologySearchForm(FlaskForm):
    protein_name = StringField('Protein name', validators=[DataRequired()])
    protein_seq = TextAreaField('Protein sequence', validators=[DataRequired(), is_protein_sequence])
    trim = BooleanField('Trim', default=False)
    ancestral = BooleanField('Ancestral', default=False)
    submit = SubmitField('Submit')

class PhyloTreeForm(FlaskForm):
    job_name = StringField('Job name', validators=[DataRequired()])
    # aln = TextAreaField('Fasta', validators=[DataRequired()])
    selection = RadioField(
        'File options', choices=[('alignment', 'sequence alignment (recommended)'),('sequence','sequences')], validators=[DataRequired()]
    )
    file = MultipleFileField('Alignment',validators=[DataRequired()])
    ancestral = BooleanField('Ancestral', default=False)
    submit = SubmitField('Submit')

class HomoModelForm(FlaskForm):
    # example_check = FormField(ExampleForm)
    job_name = StringField('Job name', validators=[DataRequired()])
    fasta = TextAreaField('Fasta', validators=[DataRequired()])

    #option 1: auto multiple template(3 tempates are used) modelling
    check_auto_multiple = BooleanField('use multiple template', default=False)

    # option 2: use user's template structure
    check_template = BooleanField('upload your own template', default=False)
    template_single = FileField('your template structure', validators=([RequiredIf(check_template=True)]))

    # option3: Homodimer modelling
    check_homodimer = BooleanField('Homodimer modelling', default=False)
    dimeT = FileField('your dimer template structure',validators=([RequiredIf(check_homodimer=True)]))
    singleT = FileField('your dimer template structure',validators=([RequiredIf(check_homodimer=True)]))
    start_residue = StringField('starting residue number', validators=([RequiredIf(check_homodimer=True)]))
    end_residue = StringField('ending residue number', validators=([RequiredIf(check_homodimer=True)]))

    submit = SubmitField('Submit')

class DockingForm(FlaskForm):
    job_name = StringField('Job name', validators=[DataRequired()])
    #fasta = TextAreaField('Fasta', validators=[DataRequired()])
    # proteins = MultipleFileField('proteins', validators=[
    #     FileRequired(), FileAllowed('pdb', 'pdb format only!')])
    # ligands = MultipleFileField('ligands', validators=[
    #     FileRequired(), FileAllowed('pdb', 'pdb format only!')])
    # ref_ligand = FileField('ref_ligand', validators=[
    #     FileRequired(), FileAllowed('pdb', 'pdb format only!')])
    proteins = MultipleFileField('proteins',validators=[DataRequired()])
    ligands = MultipleFileField('ligands', validators=[DataRequired()])
    add_h = BooleanField('add hydrogens')
    gen_3D = BooleanField('generate 3D structure')
    ref_ligand = FileField('ref_ligand')
    cof = FileField('cof')
    Rg = StringField('docking box size')
    submit = SubmitField('Submit')

class MDForm(FlaskForm):
    job_name = StringField('Job name', validators=[DataRequired()])
    job_detail = MultipleFileField('Details for the uploaded files', validators=[DataRequired()])
    proteins = MultipleFileField('Proteins', validators=[DataRequired()])
    ligands = MultipleFileField('Ligands')
    submit = SubmitField('Submit')


#############workflow_tools forms
class ScreenPanelDesForm(FlaskForm):
    job_name = StringField('Job name', validators=[DataRequired()])
    # protein_name = StringField('Protein name', validators=[DataRequired()])
    check_screened = BooleanField('Screened database')
    check_public = BooleanField('Public database')
    check_inhouse = BooleanField('In-house panel optimization', default=False)
    # protein_seq = TextAreaField('Protein sequence', validators=[RequiredIf(check_denovo=True), is_protein_sequence])
    ligands = MultipleFileField('ligands', validators=[DataRequired()])
    panel_sequences = FileField('panel_sequences', validators=([RequiredIf(check_inhouse=True)]))
    # ref_complex = FileField('ref_complex', validators=[FileRequired()])
    N_index = StringField('imine N atom index', validators=[DataRequired()])
    # cof_name = StringField('Cofactor name')
    # lig_name = StringField('Ligand name', validators=[DataRequired()])
    # check_homodimer = BooleanField('Homodimer modelling', default=False)
    # dimeT = FileField('your dimer template structure', validators=([RequiredIf(check_homodimer=True)]))
    # singleT = FileField('your dimer template structure', validators=([RequiredIf(check_homodimer=True)]))
    # start_residue = StringField('starting residue number', validators=([RequiredIf(check_homodimer=True)]))
    # end_residue = StringField('ending residue number', validators=([RequiredIf(check_homodimer=True)]))
    submit = SubmitField('Submit')

class SignatureFinderForm(FlaskForm):
    job_name = StringField('Job name', validators=[DataRequired()])
    panel_sequences = FileField('panel_sequences', validators=[DataRequired()])
    complex = FileField('ref_complex', validators=[DataRequired()])
    Rg = StringField('Rg', default='1.0')
    RMSD = StringField('RMSD', default='2.0')
    submit = SubmitField('Submit')
