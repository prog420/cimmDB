import pickle
from flask import Blueprint, render_template, request, flash
from rdkit.Chem import rdChemReactions
from application import db
from application.models import Reactions

add_bp = Blueprint(
    'add_bp', __name__,
    template_folder='templates',
    static_folder='static',
    static_url_path='/add/static'
)


@add_bp.route("/add")
def add():
    if request.args:
        print(request.args)

        values = {column: value for column, value in request.args.items()}
        reaction = rdChemReactions.ReactionFromSmarts(values['Smiles'])
        values['fingerprint'] = rdChemReactions.CreateStructuralFingerprintForReaction(reaction).ToBitString()
        new_substance = Reactions(binary_data=pickle.dumps(values))
        db.session.add(new_substance)
        db.session.commit()
        last_id_query = db.engine.execute("SELECT * FROM reactions WHERE id IN (SELECT MAX(id) FROM reactions)")
        last_id = last_id_query.fetchone()[0]
        flash(f'Record № {last_id} created')

    query = db.engine.execute("SELECT * FROM reactions")
    print(query.fetchone())
    form_fields = list(pickle.loads(query.fetchone()[1]).keys())
    form_fields.remove("Smiles")
    form_fields.remove("fingerprint")
    return render_template('add.html', form_fields=form_fields)