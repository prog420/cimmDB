const radio_buttons = document.getElementsByName("database");
const mol_dropdown = document.getElementsByName("mols_dropdown");
const reactions_dropdown = document.getElementsByName("reactions_dropdown");
const dropdowns = [mol_dropdown, reactions_dropdown];

const jsme_app = document.getElementById("jsme");

const logo_field = document.getElementById("logo_field");
const menu_elements = document.getElementById("menu_container");
const jsme_button = document.getElementById("jsme_button");

const sim_slider = document.getElementById("sim_slider")


window.addEventListener('resize', changeJSMEWidth);

sim_slider.addEventListener('input', function() {
    document.getElementById('sim_output').value = sim_slider.value;
})

for (let i=0; i < radio_buttons.length; i++) {
    // Adding events for radio buttons choosing database
    radio_buttons[i].onclick = function() {
        for (let j=0; j < dropdowns.length; j++) {
            if (dropdowns[j][0].id != `select_${this.value}`) {
                dropdowns[j][0].setAttribute('disabled', true);
            } else {
                dropdowns[j][0].removeAttribute('disabled');
            }
        }
    }
    // Event Listeners to show properties chosen from dropdowns
    dropdowns[i][0].addEventListener('change', function() {
        if (this.value == "Similarity" | this.value == "Reaction") {
            document.getElementById('sim_row').style.display = "table-cell";
            document.getElementById('conditions_row').style.display = "none";
        } else if (this.value == "Conditions") {
            document.getElementById('sim_row').style.display = "none";
            document.getElementById('conditions_row').style.display = "table-row";
        } else {
            document.getElementById('sim_row').style.display = "none";
            document.getElementById('conditions_row').style.display = "none";
        }
    });
};

function onSearch() {
    let smiles_value = document.getElementById("log").value;
    console.log(smiles_value);
    $.get(
        '/search',
        { smiles: smiles_value },
        function(data) {
            console.log(data)
            console.log('done!');
            //window.location.href = data;
        }
    )
};

function showApp() {
    jsme_app.style.display = "flex";
    jsme_button.style.display = "none";
    logo_field.style.display = "none";
    menu_elements.style.display = "none";
};

function closeApp() {
    jsme_app.style.display = "none";
    jsme_button.style.display = "";
    logo_field.style.display = "";
    menu_elements.style.display = "";
};

function changeJSMEWidth() {
    let jsme_container = document.getElementById("jsme_container").firstChild;
    jsme_container.style.width = "53vw";
    jsme_container.style.height = "50vh";
};