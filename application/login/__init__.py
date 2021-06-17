from flask import Blueprint, render_template, redirect, url_for, request, flash
from werkzeug.security import generate_password_hash, check_password_hash
from flask_login import login_user, logout_user, login_required
from flask import current_app as app
from application.models import User

# Blueprint Configuration
auth_bp = Blueprint(
    'auth_bp', __name__,
    template_folder='templates',
    static_folder='static',
    static_url_path='/auth/static'
)


@auth_bp.route("/auth", methods=['GET'])
def auth():
    """
    Homepage
    """
    return render_template("auth.html")