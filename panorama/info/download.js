function table_to_csv(source) {
    const columns = Object.keys(source.data)
    const nrows = source.get_length()
    const lines = [columns.join('\t')]

    for (let i = 0; i < nrows; i++) {
        let row = [];
        for (let j = 0; j < columns.length; j++) {// Begin at one to remove index column
            const column = columns[j]
            row.push(source.data[column][i].toString())
        }
        lines.push(row.join('\t'))
    }
    return lines.join('\n').concat('\n')
}


const filetext = table_to_csv(source)
const blob = new Blob([filetext], { type: 'text/tsv;charset=utf-8;' })

//addresses IE
if (navigator.msSaveBlob) {
    navigator.msSaveBlob(blob, filename)
} else {
    const link = document.createElement('a')
    link.href = URL.createObjectURL(blob)
    link.download = filename
    link.target = '_blank'
    link.style.visibility = 'hidden'
    link.dispatchEvent(new MouseEvent('click'))
}
