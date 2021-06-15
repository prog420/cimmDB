from werkzeug.security import generate_password_hash,  check_password_hash
from flask_login import UserMixin
from application import db


class User(UserMixin, db.Model):
    __tablename__ = 'users'

    id = db.Column(db.Integer, primary_key=True)  # primary keys are required by SQLAlchemy
    email = db.Column(db.String(100), unique=True)
    password = db.Column(db.String(100))
    name = db.Column(db.String(100), unique=True)


class Substance(db.Model):
    __tablename__ = "substances"

    id = db.Column(db.Integer, primary_key=True)
    chembl_id = db.Column(db.String(20), unique=True)
    name = db.Column(db.String(200))
    synonyms = db.Column(db.String(400))
    mol_weight = db.Column(db.REAL)
    mol_species = db.Column(db.String(40))
    mol_formula = db.Column(db.String(200))
    smiles = db.Column(db.String(2000))

