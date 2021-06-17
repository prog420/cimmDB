from flask import Blueprint, render_template
from flask_login import login_required, current_user
from flask import current_app as app

# Blueprint Configuration
main_bp = Blueprint(
    'main_bp', __name__,
    template_folder='templates',
    static_folder='static'
)


@main_bp.route("/", methods=['GET'])
def index():
    """
    Homepage
    """
    return render_template('index.html')


@main_bp.route('/profile')
@login_required
def profile():
    return render_template('profile.html', name=current_user.name)
