import mongoengine as db
from enzyme_evolver.mongo.user_models import User

class Job(db.Document):
    folder_id = db.StringField(unique=True)
    name = db.StringField()
    owner = db.ReferenceField(User)
    status = db.StringField(default='In Queue')
    type = db.StringField()
    files = db.ListField(db.StringField())
    notes = db.StringField()

    def update_status(self, new_status, print_log=False):
        self.status = new_status
        self.save()

        if print_log is True:
            print(new_status)

    def add_file(self, new_file):
        self.files.append(new_file)
        self.save()

    def update_notes(self, new_notes):
        self.notes = new_notes
        self.save()
