from flask import Blueprint, render_template, request

add_bp = Blueprint(
    'add_bp', __name__,
    template_folder='templates',
    static_folder='static',
    static_url_path='/add/static'
)


@add_bp.route("/add")
def add():
    return render_template('add.html')