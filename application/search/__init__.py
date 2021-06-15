import base64
from math import ceil
from flask import Blueprint, render_template, request
from flask import current_app as app
from CGRtools import smiles
from application import db
from application.models import Substance


# Blueprint Configuration
search_bp = Blueprint(
    'search_bp', __name__,
    template_folder='templates',
    static_folder='static',
    static_url_path='/search/static'
)


def get_description(mol, orm=False):
    """
    Construct data to send response to client
    :param mol: RowProxy object
    :return: dict
    """
    """
    smi = smiles(mol['smiles'])
    smi.clean2d()
    svg = smi.depict()
    svg_bytes = svg.encode()
    svg_base64 = base64.b64encode(svg_bytes)
    svg_base64 = svg_base64.decode("utf-8")
    """
    if orm:
        data = {'id': mol.id,
                'chembl_id': mol.chembl_id,
                'name': mol.name,
                'mol_weight': mol.mol_weight,
                'mol_species': mol.mol_species,
                'mol_formula': mol.mol_formula,
                # 'depiction': f"data:image/svg+xml;base64,{svg_base64}",
                }
    else:
        data = {'id': mol['id'],
                'chembl_id': mol['chembl_id'],
                'name': mol['name'],
                'mol_weight': mol['mol_weight'],
                'mol_species': mol['mol_species'],
                'mol_formula': mol['mol_formula'],
                # 'depiction': f"data:image/svg+xml;base64,{svg_base64}",
                }
    return data


def get_paged_results(smiles_query, mol_weight, page):
    rows_on_page = 20

    query = f"SELECT * from substances WHERE instr(smiles, '{smiles_query}') > 0 AND mol_weight > {mol_weight}"
    count_query = f'SELECT COUNT(*) FROM ({query})'
    result = db.engine.execute(query)
    total_results = db.engine.execute(count_query).fetchall()[0][0]
    num_of_pages = ceil(total_results / rows_on_page)
    for i in range(page):
        paged_res = result.fetchmany(20)
    paged_res = list(map(dict, paged_res))
    paged_res = [get_description(mol) for mol in paged_res]
    return paged_res, total_results, num_of_pages


@search_bp.route("/search", methods=['GET', 'POST'])
@search_bp.route("/search/<int:page_num>", methods=['GET', 'POST'])
def search(page_num=1):
    """
    Page with connection to database
    """

    if request.method == 'GET':
        return render_template("search.html")
    elif request.method == 'POST':
        smiles_query = request.form['smiles']
        mol_weight = int(request.form['mol_weight'])
        current_page = int(request.form['current_page'])

        #  ------------- QUERY using engine with pure SQL -------------- #

        array, total_results, num_of_pages = get_paged_results(smiles_query, mol_weight, current_page)

        #  ------------- QUERY using ORM query system + pagination -------------- #

        # Query for exact match (case INSENSITIVE)
        # query = Substance.query.filter(db.and_(Substance.smiles.like(smiles_query), Substance.mol_weight > 0))

        # Query for substring search (case INSENSITIVE)
        # query = Substance.query.filter(db.and_(Substance.smiles.contains(smiles_query),
        #                                        Substance.mol_weight > mol_weight))
        #
        # page = query.paginate(page_num, 10, False)
        # num_of_pages = page.pages
        # total_results = page.total
        # array = [get_description(mol, orm=True) for mol in page.items]

        return {'result_array': array, 'number': total_results, 'num_of_pages': num_of_pages}
