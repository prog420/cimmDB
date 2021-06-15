from flask import Blueprint, render_template
from flask import current_app as app

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