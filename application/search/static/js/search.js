const attributes = [
    'id', 
    'chembl_id',
    'mol_weight', 
    'mol_species', 
    'mol_formula',
    'depiction', 
    'smiles'];

const rows_on_page = 20;
var page = 1;

function onSearch(step) {
    let smiles_value = document.getElementById("log").value;
    let mol_weight_value = document.getElementById("mol_weight_value").textContent;
    let results = document.getElementById("result");
    let show_page = document.getElementById("page");
    let next_page = document.getElementById("next_page");
    let prev_page = document.getElementById("prev_page");
    if (step == 0){
        page = 1;
    } else {
        page += step;
    };
    if (page > 1) {
        prev_page.disabled = false;
    } else {
        prev_page.disabled = true;
    };

    $.post(
        '/search',
        { smiles: smiles_value,
          mol_weight: mol_weight_value,
          current_page: page },
        function(data) {
            let num_of_rows = data.number;
            let result_array = data.result_array;
            let num_of_pages = data.num_of_pages;
            results.innerText = num_of_rows;
            show_page.innerText = `${page}/${data.num_of_pages}`;
            if (num_of_pages > 1) {
                next_page.disabled = false;
            } else {
                next_page.disabled = true;
            };
            if (num_of_pages <= page) {
                next_page.disabled = true;
            };
            onFillTable(result_array);
        }
    )
}

function onExtendTable() {
    // После body.onload добавить необходимые ячейки в таблицу результатов для
    // последующего заполнения даннымии
    let main_table = document.getElementById('main_table')
    for (let i = 0; i < rows_on_page; i++) {
        let new_row = main_table.insertRow(-1);
        for (let j = 0; j < 5; j++) {
            let new_cell = new_row.insertCell(-1);
            let new_text = document.createTextNode("");
            new_cell.id = `${attributes[j]}_${i}`
            new_cell.appendChild(new_text);
        }
        new_row.style.display = 'none';
    }
}

function onFillTable(array) {
    // Заполнить таблицу результатов данными из БД
    for (let i = 0; i < array.length; i++) {
        let row = document.getElementById(`id_${i}`).parentElement;
        for (let j = 0; j < 5; j++) {
            let cell = document.getElementById(`${attributes[j]}_${i}`);
            cell.innerText = array[i][attributes[j]];
            /*
            if (j == 6) {
                cell.innerHTML = `<img src='${array[i][attributes[j]]}' alt='Depiction'>`;
                console.log(cell);
                cell.firstChild.style.width = '100px';
                cell.firstChild.style.height = '100px';

                //cell.firstChild.setAttribute("height", "100px");
                //cell.firstChild.setAttribute("width", "100px");
            } else {
                cell.innerHTML = array[i][attributes[j]];
            }
            */
        }
        row.style.display = '';
    }
    for (let i = array.length; i < rows_on_page; i++) {
        let row = document.getElementById(`id_${i}`).parentElement;
        row.style.display = 'none';
    }
}