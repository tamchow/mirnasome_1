const plotDivID = '{plot_id}'
const plotDiv = document.getElementById(plotDivID)
const annotationHistory = new Set()
plotDiv.on('plotly_click', data => {
    for (const point of data.points) {
        const annotateText = point.hovertext;
        const annotations = plotDiv.layout.annotations || [];
        if (annotationHistory.has(annotateText)) {
            console.log(annotationHistory)
            console.log(annotations)
            for (let i = 0; i < annotations.length; ++i) {
                if (annotations[i].text === annotateText) {
                    annotations.splice(i, 1)
                }
            }
            annotationHistory.delete(annotateText);
        } else {
            const annotation = {
                text: annotateText,
                x: point.xaxis.d2l(point.x),
                y: point.yaxis.d2l(point.y),
                opacity: 0.8,
                bgcolor: 'white',
                bordercolor: point.fullData.marker.color,
                borderwidth: 0.5,
                arrowhead: 2,
                arrowcolor: point.fullData.marker.color,
                arrowsize: 1,
                arrowwidth: 1,
                showarrow: true
            }

            annotations.push(annotation);
            annotationHistory.add(annotateText)
        }
        Plotly.relayout(plotDivID, { annotations: annotations })
    }
}).on('plotly_clickannotation', function (event, data) {
    Plotly.relayout(plotDivID, 'annotations[' + data.index + ']', 'remove');
});