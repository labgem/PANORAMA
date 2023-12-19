function filerDataBool(source, radio_buttons) {

    let selected_rows = [];
    for (let i = 0; i < source.get_length(); i++) {
        selected_rows.push(i);
    }
    for (let j = 0; j < radio_buttons.length; j++) {
        const radio_button_group = radio_buttons[j];
        const column_name = radio_button_group.name;
        const selected_value = radio_button_group.labels[radio_button_group.active];

        if (selected_value !== 'all') {
            for (let k = 0; k < selected_rows.length; k++) {
                const x = (source.data[column_name][k] === 1) === (selected_value === 'true')
                if (!x){
                    delete selected_rows[k]
                }
            }
        }
    }
    console.log(selected_rows)
    selected_rows = selected_rows.filter( function(x) {
        return x !== undefined;
    });
    console.log(selected_rows)
    const new_data = {};
    const columns = Object.keys(source.data)
    for (let l = 1; l < columns.length; l++) {// Begin at one to remove index column
        const column = columns[l]
        new_data[column] = [];
        for (let m = 0; m < selected_rows.length; m++) {
            if (l === 1 || l === columns.length - 1){
                new_data[column].push(source.data[column][selected_rows[m]]);
            } else {
                new_data[column].push(source.data[column][selected_rows[m]] === 1);
            }
        }
    }
    source.data = new_data;
    source.change.emit();
}
source.data = original_source;
source.change.emit();

filerDataBool(source, radio_buttons);