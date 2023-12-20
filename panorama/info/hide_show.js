function hide_show_columns(source, columns, checkbox_groups){
    for (let l = 1; l < columns.length; l++) {
        columns[l].visible = !!checkbox_groups.active.includes(l);
    }
    source.change.emit();
}

hide_show_columns(source, columns, checkbox_group)