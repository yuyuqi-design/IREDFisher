from flask import Flask
from redis import Redis
import rq
from enzyme_evolver.config import Config
from flask_talisman import Talisman
from flask_wtf.csrf import CSRFProtect
from flask_jsglue import JSGlue
from flask_limiter import Limiter
from flask_limiter.util import get_remote_address
from flask_security import Security, MongoEngineUserDatastore, hash_password, current_user
from flask_mail import Mail
from flask_admin import Admin
from flask_session import Session
from flask_mongoengine import MongoEngine
import datetime
from pathlib import Path


# define addons here
csrf = CSRFProtect()
jsglue = JSGlue()
limiter = Limiter(key_func=get_remote_address, default_limits=["100/minute", "1800/hour"])
talisman = Talisman(content_security_policy=False)
session = Session()
mail = Mail()
admin_ext = Admin()
db = MongoEngine()

# mongo stuff here
from enzyme_evolver.mongo.user_models import User, Role
from enzyme_evolver.mongo.workflow_models import Job
from enzyme_evolver.mongo.user_model_forms import ExtendedConfirmRegisterForm, ExtendedRegisterForm
from enzyme_evolver.app.admin import MyAdminIndexView, MyModelView
user_datastore = MongoEngineUserDatastore(db, User, Role)

# import blueprints here
from enzyme_evolver.app import main_site, workflow

def create_app(config_class=Config, use_talisman=True):
    print("Create app..")

    print(f"Path is:  {Path(__file__).parents[1]}")
    app = Flask(__name__)
    app.config.from_object(config_class)

    print("Init task queues...")
    app.redis = Redis.from_url(app.config['REDIS_URL'])
    app.task_queue = rq.Queue('tasks', connection=app.redis, default_timeout=360000*10)

    print("Init addons...")
    if use_talisman == True:
        talisman.init_app(app, content_security_policy=False)

    csrf.init_app(app)
    jsglue.init_app(app)
    session.init_app(app)
    limiter.init_app(app)
    db.init_app(app)
    mail.init_app(app)
    Security(app, user_datastore,
             confirm_register_form=ExtendedConfirmRegisterForm,
             register_form=ExtendedRegisterForm)

    print("Prepare admin views..")
    admin_ext.init_app(app, index_view=MyAdminIndexView())
    admin_ext.add_view(MyModelView(User))
    admin_ext.add_view(MyModelView(Role))
    admin_ext.add_view(MyModelView(Job))

    # Create a user to test with
    @app.before_first_request
    def create_user():
        admin = user_datastore.find_or_create_role('admin', description='admin role')

        if not user_datastore.get_user(app.config['ADMIN_EMAIL']):
            user = user_datastore.create_user(email=app.config['ADMIN_EMAIL'],
                                              password=hash_password(app.config['ADMIN_PASSWORD']),
                                              confirmed_at=datetime.datetime.now())
            user_datastore.add_role_to_user(user, admin)

    @app.context_processor
    def inject_jinja_arguments():
        inject_dict = {}
        inject_dict['login_mode'] = app.config['USE_EMAIL_CONFIRMATION']

        return inject_dict

    print("Register blueprints...")
    with app.app_context():
        app.register_blueprint(main_site.bp)
        app.register_blueprint(workflow.bp)

        return app









