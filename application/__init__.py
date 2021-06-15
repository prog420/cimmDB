from flask import Flask
from flask_sqlalchemy import SQLAlchemy

# Globally accessible libraries
db = SQLAlchemy()


def init_app():
    """
    Initialize the core application
    """
    app = Flask(__name__, instance_relative_config=False)
    app.config.from_object('config.Config')

    db.init_app(app)

    with app.app_context():
        # Import parts of app
        from . import routes
        from application.auth import auth_bp
        from application.home import home_bp
        from application.search import search_bp

        # Register Blueprints
        app.register_blueprint(auth_bp)
        app.register_blueprint(search_bp)

        return app
