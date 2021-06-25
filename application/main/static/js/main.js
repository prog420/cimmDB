const radio_buttons = document.getElementsByName("rsvp")
const mol_dropdown = document.getElementsByName("mols_dropdown")
const reactions_dropdown = document.getElementsByName("reactions_dropdown")
const dropdowns = [mol_dropdown, reactions_dropdown]
console.log(dropdowns)
const jsme_app = document.getElementById("draggable");

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
        if (this.value == "Similarity") {
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
}

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
}

function showApp() {
    jsme_app.style.display = "flex";
    $( function() {
        $( "#draggable" ).draggable();
      });
};

function closeApp() {
    jsme_app.style.display = "none";
};