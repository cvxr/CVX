// ===================================================================
// Author: Matt Kruse <matt@mattkruse.com>
// WWW: http://www.mattkruse.com/
//
// NOTICE: You may use this code for any purpose, commercial or
// private, without any further permission from the author. You may
// remove this notice from your final code if you wish, however it is
// appreciated by the author if at least my web site address is kept.
//
// You may *NOT* re-distribute this code in any way except through its
// use. That means, you can include it in your product, or your web
// site, or any other form where the code is actually being used. You
// may not put the plain javascript up on your site for download or
// include it in your javascript libraries for download. 
// If you wish to share this code with others, please just point them
// to the URL instead.
// Please DO NOT link directly to my .js files from your site. Copy
// the files to your server and use them there. Thank you.
// ===================================================================

// HISTORY
// ------------------------------------------------------------------
// December 9, 2003: Added script to the Javascript Toolbox
// December 10, 2003: Added the preProcessTrees variable to allow user
//      to turn off automatic conversion of UL's onLoad
// March 1, 2004: Changed it so if a <li> has a class already attached
//      to it, that class won't be erased when initialized. This allows
//      you to set the state of the tree when painting the page simply
//      by setting some <li>'s class name as being "liOpen" (see example)
// December 19, 2005: Added cookie persistence, and removed some things
//      that frankly I wasn't using (mostly as an exercise in understanding
//      how this worked). (Michael C. Grant)
/*
This code is inspired by and extended from Stuart Langridge's aqlist code:
        http://www.kryogenix.org/code/browser/aqlists/
        Stuart Langridge, November 2002
        sil@kryogenix.org
        Inspired by Aaron's labels.js (http://youngpup.net/demos/labels/) 
        and Dave Lindquist's menuDropDown.js (http://www.gazingus.org/dhtml/?id=109)
*/

// Automatically attach a listener to the window onload, to convert the trees
addEvent( window, "load", convertTrees );
addEvent( window, "unload", saveTrees );

function addEvent( o, e, f ) {
    if ( o.addEventListener ){ 
        o.addEventListener( e, f, true ); 
        return true; 
    } else if ( o.attachEvent ) { 
        return o.attachEvent( "on" + e, f ); 
    } else {
        return false;
    }
}

// utility function to set a global variable if it is not already set
function setDefault( name, val ) {
    if ( typeof(window[name]) == "undefined" || window[name] == null ) {
        window[name]=val;
    }
}

// Full expands a tree with a given ID
function expandTree( treeId ) {
    var ul = document.getElementById( treeId );
    if ( ul == null ) { return false; }
    expandCollapseList( ul, nodeOpenClass );
}

// Fully collapses a tree with a given ID
function collapseTree( treeId ) {
    var ul = document.getElementById( treeId );
    if ( ul == null ) { return false; }
    expandCollapseList( ul, nodeClosedClass );
}

// Performs 2 functions:
// a) Expand all nodes
// b) Collapse all nodes
function expandCollapseList( ul, cName ) {
    if ( !ul.childNodes || ul.childNodes.length == 0 ) { return false; }
    // Iterate LIs
    for ( var itemi = 0 ; itemi < ul.childNodes.length ; itemi++ ) {
        var item = ul.childNodes[itemi];
        if ( item.nodeName == "LI" ) {
            // Iterate things in this LI
            var subLists = false;
            for ( var sitemi = 0; sitemi < item.childNodes.length ; sitemi++ ) {
                var sitem = item.childNodes[sitemi];
                if ( sitem.nodeName == "UL" ) {
                    subLists = true;
                    var ret = expandCollapseList( sitem, cName );
                }
            }
            if ( subLists ) {
                item.className = cName;
            }
        } else if ( item.nodeName == "UL" ) {
        	var ret = expandCollapseList( item, cName );
       	}
    }
}

// Search the document for UL elements with the correct CLASS name, then process them
function convertTrees() {
    setDefault("treeClass","mktree");
    setDefault("nodeClosedClass","liClosed");
    setDefault("nodeOpenClass","liOpen");
    setDefault("nodeBulletClass","liBullet");
    setDefault("nodeLinkClass","bullet");
    setDefault("jsOnlyClass","jsonly");
    if ( !document.createElement ) { return; } // Without createElement, we can't do anything
    var jso = document.getElementById( jsOnlyClass );
    if ( jso != null )
    	jso.className = jsOnlyClass;
    var cndx = 0;
    var cookie = get_cookie( treeClass );
    var uls = document.getElementsByTagName( "ul" );
    for ( var uli = 0 ; uli < uls.length ; uli++ ) {
        var ul = uls[uli];
        if ( ul.nodeName == "UL" && ul.className == treeClass )
            cndx = processList( ul, cookie, cndx );
        else {
           var pn = ul.parentNode;
           if ( pn.nodeName == "DIV" && pn.className == treeClass ) {
           	   ul.className = treeClass;
               cndx = processList( ul, cookie, cndx );
           }
        }
    }
}

function saveTrees() {
    var cookie = "";
    var uls = document.getElementsByTagName( "ul" );
    for ( var uli = 0 ; uli < uls.length ; uli++ ) {
        var ul = uls[uli];
        if ( ul.nodeName == "UL" && ul.className == treeClass )
            cookie += scanList( ul );
        else {
        	var pn = ul.parentNode;
        	if ( pn.nodeName == "DIV" && pn.className == treeClass )
        	    cookie += scanList( ul );
       }
    }
    set_cookie( treeClass, cookie );
}

function scanList( ul ) {
    var cookie = "";
    if ( ul.childNodes ) {
        for ( var itemi = 0 ; itemi < ul.childNodes.length ; itemi++ ) {
            var item = ul.childNodes[itemi];
            if ( item.nodeName == "LI" ) {
                var subLists = false;
                for ( var sitemi = 0; sitemi < item.childNodes.length; sitemi++ ) {
                    var sitem = item.childNodes[sitemi];
                    if ( sitem.nodeName == "UL" ) {
                        subLists = true;
                        cookie += scanList( sitem );
                    }
                }
                if ( subLists )
                    cookie += item.className == nodeOpenClass ? "-" : "+";
            }
        }
    }
    return cookie;
}

// Process a UL tag and all its children, to convert to a tree
function processList( ul, cstr, cndx ) {
    if ( !ul.childNodes || ul.childNodes.length == 0 ) { return; }
    // Iterate LIs
    for ( var itemi = 0 ; itemi < ul.childNodes.length ; itemi++ ) {
        var item = ul.childNodes[itemi];
        if ( item.nodeName == "LI" ) {
            // Iterate things in this LI
            var subLists = false;
            for ( var sitemi = 0; sitemi < item.childNodes.length; sitemi++ ) {
                var sitem = item.childNodes[sitemi];
                if ( sitem.nodeName == "UL" ) {
                    subLists = true;
                    cndx = processList( sitem, cstr, cndx );
                }
            }
            var s = document.createElement("SPAN");
            var t = '\u00A0'; // &nbsp;
            s.className = nodeLinkClass;
            if (subLists) {
                // This LI has UL's in it, so it's a +/- node
                if ( cndx < cstr.length ) {
                    thisC = cstr.substring( cndx, cndx + 1 ); 
                    item.className = thisC == "-" ? nodeOpenClass : nodeClosedClass;
                    ++cndx;
                } else if ( item.className == null || item.className == "" )
                    item.className = nodeClosedClass;
                // If it's just text, make the text work as the link also
                if ( item.firstChild.nodeName=="#text" ) {
                    t = t + item.firstChild.nodeValue;
                    item.removeChild( item.firstChild );
                }
                s.onclick = function () {
                    this.parentNode.className = ( this.parentNode.className == nodeOpenClass ) ? nodeClosedClass : nodeOpenClass;
                    return false;
                }
            } else {
                // No sublists, so it's just a bullet node
                item.className = nodeBulletClass;
                s.onclick = function () { 
                    return false; 
                }
            }
            s.appendChild( document.createTextNode( t ) );
            item.insertBefore( s, item.firstChild );
        }
    }
    return cndx;
}

function get_cookie( name ) {
    var returnvalue = "";
    cookie = document.cookie;
    if ( cookie.length > 0 ) {
        var search = name + "=";
        offset = cookie.indexOf( search );
        if ( offset != -1 ) {
            offset += search.length;
            end = cookie.indexOf( ";", offset );
            if ( end == -1 ) end = cookie.length;
            returnvalue = unescape( cookie.substring( offset, end ) );
        }
    }
    return returnvalue;
}

function set_cookie( name, value ) {
    document.cookie = name + "=" + value + ";";
}
