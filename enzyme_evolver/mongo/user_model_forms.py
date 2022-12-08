from flask_security import RegisterForm, ConfirmRegisterForm
from wtforms import StringField, PasswordField, BooleanField, SubmitField
from wtforms.validators import DataRequired, EqualTo
from flask_wtf import FlaskForm
from wtforms.fields.html5 import EmailField

class ExtendedConfirmRegisterForm(ConfirmRegisterForm):
    password_confirm = PasswordField('Retype Password', validators=[EqualTo('password', message='RETYPE_PASSWORD_MISMATCH')])
    email_opt_in = BooleanField('I am willing to be contacted (infrequently) by email',
                                default='checked')

class ExtendedRegisterForm(RegisterForm):
    password_confirm = PasswordField('Retype Password', validators=[EqualTo('password', message='RETYPE_PASSWORD_MISMATCH')])
    email_opt_in = BooleanField('I am willing to be contacted (infrequently) by email',
                                default='checked')

class UserProfileForm(FlaskForm):
    email_opt_in = BooleanField('I am willing to be contacted (infrequently) by email',
                                default='checked')
    submit = SubmitField('Update information')


