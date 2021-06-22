from flask import Blueprint, render_template
from flask_login import login_required, current_user

# Blueprint Configuration
main_bp = Blueprint(
    'main_bp', __name__,
    template_folder='templates',
    static_folder='static',
    static_url_path='/main/static'
)

conditions = {'Acid hydrolysis',
              'Autoclave',
              'Darkness',
              'Electrolysis',
              'Green chemistry',
              'Heating',
              'Heating / reflux',
              'Inert atmosphere',
              'Irradiation',
              'Reflux',
              'Schlenk technique',
              'Sealed tube',
              'UV-irradiation',
              'cooling',
              'infrared irradiation',
              'microwave irradiation',
              'sonication'}


@main_bp.route("/", methods=['GET'])
def index():
    """
    Homepage
    """
    return render_template('main.html', name=current_user, conditions=conditions)


@main_bp.route('/profile')
@login_required
def profile():
    return render_template('profile.html', name=current_user.name)
