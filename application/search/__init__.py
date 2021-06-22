import time
from math import ceil
from flask import Blueprint, render_template, request, make_response
from flask_restful import Resource
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from application import db
from application.models import Substance


# Blueprint Configuration
search_bp = Blueprint(
    'search_bp', __name__,
    template_folder='templates',
    static_folder='static',
    static_url_path='/search/static'
)


class Search(Resource):
    data = []
    similarity = 0.0

    def get(self):
        pass

    def put(self):
        self.data.append(request.form['data'])
        return self.data

    def post(self):
        """
        Function to handle user POST query
        fingerprints - array of ExplicitBitVectors
        similarity score - Tanimoto
        Estimation performs with BulkTanimotoSimilarity - fastest method for RDKit DataStructures
        Algorithm of similarity search:
            1. Obtaining available molecules
            2. Forming array of bitvectors based on bitstrings
            3. Calculating similarity for all pairs "query fingerprint - molecule fingerprint"
            4. Forming array of enumerated values to get needed molecules
            4. Filtering array with needed accuracy
            5. Obtaining molecules with required accuracy
        """
        user_query = request.form.get('smiles')
        query_type = request.form.get('rsvp')
        subquery_type = request.form.get('rsvp_dropdown')
        if subquery_type == 'Similarity':
            self.similarity = int(request.form.get('similarity_value'))/100
        query_mol = Chem.MolFromSmiles(user_query)
        query_fp = FingerprintMols.FingerprintMol(query_mol).ToBitString()
        query_fp = query_fp + '0'*(2048 - len(query_fp))
        query_fp = DataStructs.cDataStructs.CreateFromBitString(query_fp)
        query = Substance.query.all()
        if user_query == "":
            self.data = query

        else:
            start_time = time.time()
            fingerprints = [DataStructs.cDataStructs.CreateFromBitString(mol.fingerprint) for mol in query]
            similarities = list(enumerate(DataStructs.BulkTanimotoSimilarity(query_fp, fingerprints)))
            filtered_similarities = list(filter(self.check_similarity, similarities))
            query = [query[indice[0]] for indice in filtered_similarities]
            self.data = similarities
            print("--- %s seconds ---" % (time.time() - start_time))
        if query:
            molecules = [{'id': mol.id,
                          'chembl_id': mol.chembl_id,
                          'mol_weight': mol.mol_weight,
                          'hba_lipinski': mol.hba_lipinski,
                          'hbd_lipinski': mol.hbd_lipinski,
                          'mol_species': mol.mol_species,
                          'mol_formula': mol.mol_formula}
                         for mol in query[0:10]]
            amount = len(query)
        else:
            molecules = [{}]
            amount = 0
        headers = {'Content-Type': 'text/html'}
        return make_response(render_template('search.html', molecules=molecules, amount=amount), 200, headers)

    def check_similarity(self, sim_value):
        """
        Check if similarity of queried fingerprint and fingerprint from database
        is high enough to keep it
        :param sim_value: Tuple(enumerate, similarity_value)
        :return: Bool
        """
        if sim_value[1] >= self.similarity:
            return True
        else:
            return False


def get_description(mol, orm=False):
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


"""
@search_bp.route("/search", methods=['GET', 'POST'])
@search_bp.route("/search/<int:page_num>", methods=['GET', 'POST'])
def search(page_num=1):
    ""
    Page with connection to database
    ""

    smiles_query = request.form.get('smiles')
    print(smiles_query)
    query = Substance.query.all()
    print(type(query))
    print(dir(query))
    result = []
    for mol in query:
        result.append(mol.smiles)
    print(result)
    return render_template('search.html')
    # current_page = int(request.form['current_page'])


    #  ------------- QUERY using engine with pure SQL -------------- #

    # array, total_results, num_of_pages = get_paged_results(smiles_query, mol_weight, current_page)

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
"""
