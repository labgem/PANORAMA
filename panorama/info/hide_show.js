// const {ColumnDataSource} = require("@bokeh/bokehjs/lib/models/widgets/main");
/**
 * @param source ColumnDataSource
 * @param {number} index
 */
function hide_show_column(columns, index) {
    console.log(index)
    let column_to_hide_show = columns[index]
    console.log(column_to_hide_show)
    console.log(column_to_hide_show.visible)
    column_to_hide_show.visible = !column_to_hide_show.visible;
}

hide_show_column(columns, index)