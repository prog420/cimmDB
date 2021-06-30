import sqlite3
from datetime import timedelta
from flask import Flask, g
from flask_sqlalchemy import SQLAlchemy
from flask_restful import Api
from flask_login import LoginManager
from flask_session import Session

# Globally accessible libraries
db = SQLAlchemy()
sess = Session()

login_manager = LoginManager()
login_manager.login_view = 'auth_bp.login.html'


def init_app():
    """
    Initialize the core application
    """
    app = Flask(__name__, instance_relative_config=False)
    app.config.from_object('config.Config')

    app.permanent_session_lifetime = timedelta(minutes=5)

    db.init_app(app)
    login_manager.init_app(app)
    sess.init_app(app)

    with app.app_context():
        # Import parts of app
        # from . import routes
        from application.login import auth_bp
        from application.main import main_bp
        from application.search import search_bp, Search
        from application.add import add_bp
        from application.models import User

        @login_manager.user_loader
        def load_user(user_id):
            # Since the user_id is just the primary key of our user table
            # use it in the query for the user
            return User.query.get(int(user_id))

        search_routes = [
            "/search",
            "/search/<int:page_num>",
        ]

        api = Api(search_bp)
        # api.add_resource(Search, "/search/<int:page_num>", endpoint='search')
        api.add_resource(Search, "/search/<int:page_num>", *search_routes)

        # Register Blueprints
        app.register_blueprint(main_bp)
        app.register_blueprint(auth_bp)
        app.register_blueprint(search_bp)
        app.register_blueprint(add_bp)

        return app

