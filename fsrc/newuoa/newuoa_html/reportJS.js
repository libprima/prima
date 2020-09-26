$(window).scroll(function () {
    var filename = (document.location.pathname.match(/[^\/]+$/)[0]).split('.')[0]   //Gets the filename without .html extension (ie index.html == index)

    if ($("#nav-dropdown-content").css('display') == 'block'){  //Stops sticky nav collapse if Contents is displayed (prevents everything jumping around)
        if ($(document).scrollTop() > 40) {
            sessionStorage.setItem('navbarScrollHold', 'expanded');    //Stores variable to tell page navbar is expanded (but should be collapsed)
        }
        else if ($(document).scrollTop() <= 40) {
            sessionStorage.setItem('navbarScrollHold', 'collapsed');    //Stores variable to tell page navbar is collapsed (but should be expanded)
        }
        return;
    }
    else if ($(document).scrollTop() <= 40) {
        if ($('#main-nav').hasClass('title-scrolled-nav')) {
            if ($(document).scrollTop() <= 3) {
                $('#main-nav').removeClass('title-scrolled-nav');     // Changes stylings based on scrolling, special version for user generated title page
                $('#main-nav').css('height', '113px');
            //$('#titlepg-main-content').css('margin-top', '95px');
            window.scrollTo(0,0);
            }
        } else {
            var scrollSwitch = 0;
            if ($('#main-nav').hasClass('scrolled-nav')) {
                scrollSwitch = 1;
            }
            $('#main-nav').removeClass('scrolled-nav');     // Changes different stylings based on user scrolling up to top 40px
            $('#main-content').removeClass('scrolled-nav');
            if (scrollSwitch == 1) {
                window.scrollTo(0,0);
            }
        }
    } 
    else {
        if (filename == "titlepage") {          // Specific styling for user generated title page
            var titleScrollSwitch = 0;
            if (!$('#main-nav').hasClass('title-scrolled-nav')) {
                titleScrollSwitch = 1;
            }
            $('#main-nav').addClass('title-scrolled-nav');
            $('#main-nav').css('height', '0px');
            //$('#titlepg-main-content').css('margin-top', '41px');       // Snaps page down to title page when navbar vanishes on scroll
            if (titleScrollSwitch == 1) {
                window.scrollTo(0,4);
                setTimeout( function() {
                    window.scrollTo(0,4);
                }, 200);
            }
        } else {
            $('#main-nav').addClass('scrolled-nav');        // Changes different stylings based on user scrolling down past 40px
            $('#main-content').addClass('scrolled-nav');
        }
    } 
}).one('scroll', function() {
    $('#dropdown-buffer').css("display", "none");
    $('.dropdown-content').css("display","none");   // Hides initially displayed dropdown if user scrolls page at all
    }); 

// Function checks every .2s to see if navbar is being held and can now be collapsed/expanded (respectively)
setInterval(function() {
    var filename = (document.location.pathname.match(/[^\/]+$/)[0]).split('.')[0]   //Gets the filename without .html extension (ie index.html == index)

    if (sessionStorage.getItem('navbarScrollHold') == 'expanded' && $("#nav-dropdown-content").css('display') == 'none') {
        if ($(document).scrollTop() > 40) {
            if (filename == "titlepage") {
                $('#main-nav').addClass('title-scrolled-nav');
                //$('#titlepg-main-content').css('margin-top', '41px');
            } else {
                $('#main-nav').addClass('scrolled-nav');        // Changes different stylings based on user scrolling down past 40px
            }
        }
        sessionStorage.setItem('navbarScrollHold', 'no');
    }
    else if (sessionStorage.getItem('navbarScrollHold') == 'collapsed' && $("#nav-dropdown-content").css('display') == 'none') {
        if ($(document).scrollTop() <= 40) {
            if ($('#main-nav').hasClass('title-scrolled-nav')) {
                $('#main-nav').removeClass('title-scrolled-nav');     // Changes stylings based on scrolling, special version for user generated title page
                //$('#titlepg-main-content').css('margin-top', '95px');
            } else {
                $('#main-nav').removeClass('scrolled-nav');     // Changes different stylings based on user scrolling up to top 40px
                window.scrollTo(0,0);
            }
        }
        sessionStorage.setItem('navbarScrollHold', 'no');
    }
}, 200);

// Function adds a hover buffer to the table of contents, so slightly missing the edge will not close the ToC
$(document).ready( function () {
$(".dropdown").mouseover(
    function() {
    var height = $('#nav-dropdown-content').height();
    $('#dropdown-buffer').css('height', (height+50) + 'px')
});
});

// Left and right table scroll buttons
function scrollTableRight() {
    var scroll = $('#main-content').scrollLeft();
    $('#main-content').scrollLeft(scroll + 50);
};

function scrollTableLeft() {
    var scroll = $('#main-content').scrollLeft();
    $('#main-content').scrollLeft(scroll - 50);
};

// Functions to detect whether table is wider than viewport, and if so, make buttons appear
function tableWidth() {
    var tableWidth = $('#metrics-table').width();
    var viewportWidth = $(window).width();
    // alert(tableWidth + "," + viewportWidth);
    if (tableWidth+20 > viewportWidth) {
        $('.navbar #rightspacer button').css('display', 'block');
        $('#imgRight').css('display', 'block');
        $('#imgLeft').css('display', 'block');
    } else if (tableWidth+20 < viewportWidth) {
        $('#imgRight').css('display', 'none');
        $('#imgLeft').css('display', 'none');
        $('.navbar #rightspacer button').css('display', 'none');
    }

};
$(document).ready( function() {
    tableWidth();

});
$(window).resize( function() {
    tableWidth() 
});


// Removes shadow from initial navbar (done this way so shadow still appears when JS not supported)
$(document).ready( function() {
    $('#lower-nav').css('box-shadow', 'none');
});

// Makes the dropdown which is displayed when page is first opened disapper after clicking anywhere
$(document).one('click', function() {
  $('.dropdown-content').css("display","none");
  $('#dropdown-buffer').css("display", "none");
});

// Returns the current page to the top
function pageToTop() {
    document.body.scrollTop = 0;
    document.documentElement.scrollTop = 0;
};

// Checks whether this is the first time a user is viewing the report in a given session, and if it is, displays the Table of Contents dropdown
// Note that this uses sessionStorage, which is cleared when browser window is closed. localStorage could be used to persist forever
$(document).ready(
    function() {
if (!sessionStorage.getItem('viewedUnderstandReport')){
        $('.dropdown-content').css("display", "block");
        var height = $('#nav-dropdown-content').height();
        $('#dropdown-buffer').css("display", "block");
        $('#dropdown-buffer').css('height', "100vh");
        sessionStorage.setItem('viewedUnderstandReport', 'yes');
    }
});

// NavBar: navbar() function is called from html page on load, places upperNav and table of contents (var ToC
//          in file table_of_contents.js) into html divs by template literal
function navbar() {
    $('#lower-nav').prepend(ToC);
    var htmlToc = document.getElementById("main-dropdown-html");
    htmlToc.parentNode.removeChild(htmlToc);
};


