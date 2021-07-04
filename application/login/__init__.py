import datetime
from flask import Blueprint, render_template, redirect, url_for, request, flash
from werkzeug.security import generate_password_hash, check_password_hash
from flask_login import login_user, logout_user, login_required
from flask import current_app as app
from application import db
from application.models import User


# Blueprint Configuration
auth_bp = Blueprint(
    'auth_bp', __name__,
    template_folder='templates',
    static_folder='static',
    static_url_path='/login/static'
)


@auth_bp.route("/login")
def login():
    """
    Homepage
    """
    return render_template("login.html")


@auth_bp.route('/login', methods=['POST'])
def login_post():
    email = request.form.get('email')
    password = request.form.get('password')
    remember = True if request.form.get('remember') else False

    user = User.query.filter_by(email=email).first()

    # Check if user actually exists
    # Take the user supplied password, hash it, and compare it to
    # the hashed password in database
    if not user or not check_password_hash(user.password, password):
        flash('Please check your login details and try again.')
        # If user doesn't exist or password is wrong, reload the page
        return redirect(url_for('auth_bp.login'))

    # If the above check passes, then we know the user has the right credentials
    login_user(user, remember=remember, duration=datetime.timedelta(days=31))
    return redirect(url_for('main_bp.profile'))


@auth_bp.route('/signup')
def signup():
    return render_template('signup.html')


@auth_bp.route('/signup', methods=['POST'])
def signup_post():
    email = request.form.get('email')
    name = request.form.get('name')
    password = request.form.get('password')

    # If this returns a user, then the email already exists in database
    # If a user is found, we want to redirect back to signup page so user can try again
    if User.query.filter_by(email=email).first():
        flash('Email address already exists.')
        return redirect(url_for('auth_bp.signup'))
    elif User.query.filter_by(name=name).first():
        flash('Name of user already exists.')
        return redirect(url_for('auth_bp.signup'))

    # Create new user with the form data. Hash the password so plaintext version won't be saved.
    new_user = User(email=email, name=name, password=generate_password_hash(password, method='sha256'))

    # Add the new user to the database
    db.session.add(new_user)
    db.session.commit()

    return redirect(url_for('auth_bp.login'))


@auth_bp.route('/logout')
@login_required
def logout():
    logout_user()
    return redirect(url_for('main_bp.index'))