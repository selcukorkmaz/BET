var amountScrolled = 300;

$(window).scroll(function() {
    if ( $(window).scrollTop() > amountScrolled ) {
        $('a.back-to-top').fadeIn('slow');
    } else {
        $('a.back-to-top').fadeOut('slow');
    }
});


$(document).keyup(function(event) {
    if ($("#PdbId").is(":focus") && (event.keyCode == 13)) {
        $("#startAnalysis").click();
    }
});


$(document).keyup(function(event) {
    if ($("#PDFPdbId").is(":focus") && (event.keyCode == 13)) {
        $("#startAnalysis").click();
    }
});


$(document).keyup(function(event) {
    if ($("#PDFPubID").is(":focus") && (event.keyCode == 13)) {
        $("#startAnalysis").click();
    }
});

$(document).ready(function() {
    $('#TextMining').DataTable();
} );