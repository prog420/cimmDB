"""Flask configuration."""
from redis import Redis


class Config:
    """Set Flask configuration from environment variables."""

    FLASK_APP = 'wsgi.py'
    FLASK_ENV = 'development'
    SECRET_KEY = '083aa2439117d25e0bbaf697b97c2189ded70bb97cedc170'

    # Redis
    SESSION_TYPE = 'redis'
    SESSION_REDIS = Redis(host="redis-18039.c278.us-east-1-4.ec2.cloud.redislabs.com", port="18039",
                          password="D8qYzmrPaNPwudZK1t5Mh3z0SOJqhEpt")

    TESTING = True
    DEBUG = True

    # Static Assets
    STATIC_FOLDER = 'static'
    TEMPLATES_FOLDER = 'templates'

    # Flask-SQLAlchemy
    # SQLALCHEMY_DATABASE_URI = "postgresql://postgres:password@localhost/test"
    SQLALCHEMY_DATABASE_URI = "sqlite:///db.sqlite"
    SQLALCHEMY_ECHO = False
    SQLALCHEMY_TRACK_MODIFICATIONS = False
