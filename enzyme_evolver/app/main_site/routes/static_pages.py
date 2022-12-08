from enzyme_evolver.app.main_site import bp
from flask import render_template,send_file
from pathlib import Path

@bp.route('/', methods=['GET'])
def enzyme_evolver():
    return render_template('enzyme_evolver_home.html')
@bp.route('/IREDFisher', methods=['GET'])
def fisher_home():
    return render_template('fisher_home.html')

@bp.route('/SignatureFinder', methods=['GET'])
def Sigfinder_home():
    return render_template('Sigfinder_home.html')

@bp.route('/ComputHub', methods=['GET'])
def computhub_home():
    return render_template('computhub_home.html')

@bp.route('/IREDFisher/contact', methods=['GET'])
def contact():
    return render_template('contact.html')

@bp.route('/ComputHub/contact', methods=['GET'])
def contact_computhub():
    return render_template('contact_computhub.html')

@bp.route('/enzyme_fisher/tutorials', methods=['GET'])
def tutorials():
    return render_template('tutorials.html')

@bp.route('/enzyme_fisher/tutorial_files', methods=['GET'])
def tutorial_files():
    file=f'{Path.cwd()}/enzyme_evolver/database/tutorial.zip'
    return send_file(file, attachment_filename='tutorial.zip')
@bp.route('/enzyme_fisher/mannual', methods=['GET'])
def mannual():
    file=f'{Path.cwd()}/enzyme_evolver/database/EnzymeFisher_mannual.pdf'
    return send_file(file, attachment_filename='EnzymeFisher_mannual.pdf')


@bp.route('/cookie_policy', methods=['GET'])
def cookie_policy():
    return render_template('cookie_policy.html')

@bp.route('/terms', methods=['GET'])
def terms():
    return render_template('terms.html')
