/*
    Rpad_editor.js  --  creates the "lite" editor that works in IE or Mozilla

    by Tom Short, EPRI Solutions, Inc.., tshort@eprisolutions..com

	(c) Copyright 2005. by EPRI Solutions, Inc..

	Permission is hereby granted, free of charge, to any person obtaining a
	copy of this software and associated documentation files (the "Software"),
	to deal in the Software without restriction, including without limitation
	the rights to use, copy, modify, merge, publish, distribute, sublicense,
	and/or sell copies of the Software, and to permit persons to whom the
	Software is furnished to do so, subject to the following conditions:

	The above copyright notice and this permission notice shall be included
	in all copies or substantial portions of the Software.

	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
	OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
	ITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
	AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
	LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
	FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
	DEALINGS IN THE SOFTWARE.
*/

HTMLArea._addEvent(window, "load", startRpad);  // With this, we don't have to use <body onload="startRpad"> in the HTML

var _Rpad_editor = null;
var _Rpad = null;

// This is the main initialization routine.
function initEditor() {
  add_div_wrapper();
                                                                                                                                                                                       
// set up the menu bar and status area                                                                                                                                               
                                                                                                                                                                                       
  var el = document.createElement("div"); // blank div to take up space
  el.height = "30px"; 
  el.innerHTML = "<div>&nbsp</div>";
  document.body.insertBefore(el,document.body.firstChild);                                                                                                                             
  var el = document.createElement("div"); // the top bar
  el.setAttribute('contentEditable',false);
  el.id = "Rpad_topBar";

//  el.setAttribute('unselectable','On');
  el.innerHTML = '<table height="30px" width="100%" contentEditable="false"><td><div id="Rpad_menuBar" unselectable="on" style="-moz-user-select:none;"></div></td><td align="right"><div id="Rpad_statusArea" style="fontsize:3 -moz-user-focus: normal; -moz-user-select: normal;"></div></td></table>';
  document.body.insertBefore(el,document.body.firstChild);                                                                                                                             
  document.body.setAttribute('unselectable','On');
                                                                                                                                                                                       
  cmDraw('Rpad_menuBar', Rpad_menuBar, 'hbr', cmThemeOffice, 'ThemeOffice');                                                                                                           

  _Rpad_editor = new HTMLArea();
  _Rpad_editor._doc = document;
  _Rpad_editor._statusArea = document.getElementById("Rpad_statusArea");
  _Rpad_editor._menuBar = document.getElementById("Rpad_menuBar");
  _Rpad = document.getElementById("Rpad_div");

  window.onbeforeunload = stopRpad;
//  window.beforeunload = stopRpad;
//  window.onunload = stopRpad;

//  document.write('<style type="text/css">@import url(editor/Rpad_page.css);</style>');

  HTMLArea._addEvents
//    (document, ["keydown", "keypress", "mousedown", "mouseup", "drag"],
    (_Rpad, ["keydown", "keypress", "mousedown", "mouseup", "drag"],
        function (event) {
            return _editorEvent(event);
        });
}

// Add a textarea wrapper around the whole document body.
// This is somewhat cludgy, but htmlarea is a replacement for textarea's and not a general page editor.
// Also, watch out that the wrapper must skip over the first document element, which should be the table
// that holds the menu bar.
function add_div_wrapper() {
    var d = document.body;
    var div = document.createElement('div');
    div.setAttribute('contentEditable',true);
    div.id = "Rpad_div";
    d.insertBefore(div,d.firstChild);
    while (d.firstChild.nextSibling) // move the body's children
      div.appendChild(d.firstChild.nextSibling);
}  
   
// This is a modified version of the resize_editor in htmlarea.js used for fullscreen operation.
function resize_editor() {  // resize editor to fix window
  var newHeight, newWidth;
  if (document.all) {
    // IE
    newHeight = document.body.offsetHeight - _Rpad_editor._menuBar.offsetHeight - 15;
    newWidth = document.body.offsetWidth - 6;
    if (newHeight < 0) { newHeight = 0; }
  } else {
    // Gecko
    newHeight = window.innerHeight - _Rpad_editor._menuBar.offsetHeight - 24;
    newWidth = window.innerWidth - 15;
  }
  _Rpad.style.height = newHeight + "px";
  _Rpad.style.width = newWidth + "px";
}  

function _isInRpadInput() {
	var p = _Rpad_editor.getParentElement();
    var result = false;
	while (p && (p.nodeType == 1) && (p.tagName != 'BODY')) {
        result = result || 
                    (p.className == "Rpad_input" && p.tagName != 'INPUT' && p.tagName != 'SELECT' && p.tagName != 'TEXTAREA');
		p = p.parentNode;
	}
	return result;
};


/** Event handler.
 * This function also handles key bindings. */
function _editorEvent(ev) {
	var editor = this;
	var keyEvent = (HTMLArea.is_ie && ev.type == "keydown") || (ev.type == "keypress");
	if (keyEvent && ev.ctrlKey) {
		var sel = null;
		var range = null;
		var key = String.fromCharCode(HTMLArea.is_ie ? ev.keyCode : ev.charCode).toLowerCase();
		var cmd = null;
		var value = null;
		switch (key) {
		    case 'a':
			if (!HTMLArea.is_ie) {
				// KEY select all
				sel = this._getSelection();
				sel.removeAllRanges();
				range = this._createRange();
				range.selectNodeContents(document.body);
				sel.addRange(range);
				HTMLArea._stopEvent(ev);
			}
			break;

			// simple key commands follow

		    case 'b': cmd = "bold"; break;
		    case 'i': cmd = "italic"; break;
		    case 'u': cmd = "underline"; break;
		    case 's': Rpad_save(); HTMLArea._stopEvent(ev); break;
		    case 'l': cmd = "justifyleft"; break;
		    case 'e': cmd = "justifycenter"; break;
//		    case 'r': cmd = "justifyright"; break;
		    case 'j': cmd = "justifyfull"; break;
		    case 'z': cmd = "undo"; break;
		    case 'y': cmd = "redo"; break;
//		    case 'v': cmd = "paste"; break;

		    case '0': cmd = "killword"; break;

			// headings
		    case '1':
		    case '2':
		    case '3':
		    case '4':
		    case '5':
		    case '6':
			cmd = "formatblock";
			value = "h" + key;
			if (HTMLArea.is_ie) {
				value = "<" + value + ">";
			}
			break;
		}
		if (cmd) {
			// execute simple command
			_Rpad_editor._doc.execCommand(cmd, false, value);
			HTMLArea._stopEvent(ev);
		}
	}
	else if (keyEvent) {
		// other keys here
		switch (ev.keyCode) {
		    case 13: // KEY enter
         		// if (HTMLArea.is_ie) {
                if (_isInRpadInput()) {
         		  _Rpad_editor.insertHTML("<br />");
         		  HTMLArea._stopEvent(ev);
         		}
			break;

		    case 8: // BACKSPACE key - disable for Mozilla/Mozile
                 if (!HTMLArea.is_ie && ev.originalTarget.nodeName != "INPUT" && ev.originalTarget.nodeName != "TEXTAREA") { 
				   HTMLArea._stopEvent(ev);
                 }
			break;
		    case 0x78: // F9 key
                    var isF9 = (!HTMLArea.is_ie || ev.type == "keydown") ;
                    if (isF9) { // This test distinguishes between F9 & 'x'
//                    if (ev.type == "keydown") { // This test distinguishes between F9 & 'x'
                      Rpad_calculate();
				  HTMLArea._stopEvent(ev);
                    }
    			break;
		    case 0x77: // F8 key
                    var isF8 = (!HTMLArea.is_ie || ev.type == "keydown") ;
                    if (isF8) { // This test distinguishes between F8 & 'x'
                      Rpad_calculate_selection();
				  HTMLArea._stopEvent(ev);
                    }
    			break;
		    case 0x71: // F2 key
                    var isF2 = (!HTMLArea.is_ie || ev.type == "keydown") ;
                    if (isF2) { // This test distinguishes between F2 & 'x'
                      _Rpad_editor.RpadToggleCodeVisibility();
				  HTMLArea._stopEvent(ev);
                    }
    			break;
		}
	}
//editor.updateStatusArea(ev);
	// update the toolbar state after some time
	if (_Rpad_editor._timerToolbar) {
		clearTimeout(_Rpad_editor._timerToolbar);
	}
    var target = HTMLArea.is_ie ? ev.srcElement : ev.target;
	_Rpad_editor._timerToolbar = setTimeout(function() {
		_Rpad_editor._updateStatusArea(target);
		_Rpad_editor._timerToolbar = null;
	}, 50);
};

HTMLArea.prototype._updateStatusArea = function(node) {
  var parent = this.getParentElement();
  // try an INPUT or SELECT
  if (node.tagName == 'SELECT' || node.tagName == 'INPUT' || node.tagName == 'TEXTAREA') 
    parent = node;
  window.status = parent.tagName+parent.className + parent.name;

  var str = ""
  var el = document.createElement("div");
  el.style.fontSize = "12";
  el.className = "statusArea";
  if (parent.getAttribute('Rpad_type') != null) {
    var select = document.createElement("select");
    select.style.fontSize = "12";
    var opts=new Array("R script","R string","R variable","file","server shell script","javascript");
    var vals=new Array("R","Rstring","Rvariable","file","shell","javascript");
    for (i=0;i<opts.length;i++) {
      var o = document.createElement("OPTION");
      var t = document.createTextNode(opts[i]);
      o.setAttribute("value",vals[i]);
      o.appendChild(t);
      select.appendChild(o);
      if (vals[i] == parent.getAttribute("rpad_type"))
        select.selectedIndex = i;
    }
    select.onchange = function() {parent.setAttribute("rpad_type",this.value);};
    el.appendChild(select);
  }
  var canHaveName = parent.tagName.toLowerCase() == "input" || 
                    parent.tagName.toLowerCase() == "select" || 
                    (
                      (parent.getAttribute("Rpad_type") == "file" || 
                       parent.getAttribute("Rpad_type") == "Rstring" || 
                       parent.getAttribute("Rpad_type") == "Rvariable"));
  if (canHaveName) {
    var inp = document.createElement("input");
    inp.ID = "statusarea_name";
    inp.style.fontSize = "12";
//    inp.value = parent.getAttribute("name");
    inp.value = parent.getAttribute("name");
    inp.onchange = function() {parent.setAttribute("name",this.value);};
    inp.onblur = function() {parent.setAttribute("name",this.value);};
    el.appendChild(inp);
  }
  var canHaveValue = parent.tagName.toLowerCase() == "input" &&
                     parent.type.toLowerCase() == "radio";
  if (canHaveValue) {
    var inp = document.createElement("input");
    inp.ID = "statusarea_value";
    inp.style.fontSize = "12";
    inp.value = parent.getAttribute("value");
    inp.onchange = function() {parent.setAttribute("value",this.value);};
    inp.onblur = function() {parent.setAttribute("value",this.value);};
    el.appendChild(inp);
  }
  this._statusArea.parentNode.replaceChild(el,this._statusArea)
  this._statusArea = el;
};


var Rpad_menuBar =
[
	[null, 'Rpad', null, null, null,
		['<img class="seq1" src="editor/JSCookMenu/ThemeOffice/save.gif" /><img class="seq2" src="editor/JSCookMenu/ThemeOffice/saveshadow.gif" />', 'Save', 'javascript:Rpad_save()', null, null],
		[null, 'Save As', 'javascript:Rpad_save_as()', null, null],
		_cmSplit,
   		[null, 'Toggle code visibility [F2]', "javascript:_Rpad_editor.toggleCodeVisibility()", '_self', null],
   		[null, 'Unhide all code', "javascript:RpadUnhideElements()", '_self', null],
	],
	_cmSplit,
	[null, 'Insert', null, null, null,
		[null, 'Rpad section', "javascript:Rpad_insert_Rpad('R')", '_self', null],
		[null, 'Rpad textarea', "javascript:Rpad_insert_Rpad_textarea('file')", '_self', null],
		[null, 'Rpad span', "javascript:Rpad_insert_Rpad_span('R')", '_self', null],
		[null, 'Rpad input box', "javascript:Rpad_insert_Rpad_input('Rvariable')", '_self', null],
		_cmSplit,
		[null, 'Radio button', "javascript:Rpad_insert_radio_input()", '_self', null],
		[null, 'Check box', "javascript:Rpad_insert_check_input()", '_self', null],
    ],
	_cmSplit,
	[null, 'Calculate [F9]', "javascript:Rpad_calculate()", null, null, null],
	_cmSplit,
	[null, 'Help', null, null, null,
		['<img class="seq1" src="editor/JSCookMenu/ThemeOffice/help.gif" /><img class="seq2" src="editor/JSCookMenu/ThemeOffice/helpshadow.gif" />', 'Rpad', 'BasicDocumentation.html', '_blank', null],
		[null, 'Rpad homepage', 'http://www.rpad.org', '_blank', null],
		[null, 'R homepage', 'http://www.r-project.org', '_blank', null],
		[null, 'R documentation', 'R/doc/html/index.html', '_blank', null],
		[null, 'R intro', 'R-intro.Rpad', '_blank', null],
		[null, 'R function finder', 'SearchRKeywords.Rpad', '_blank', null]
	]
];

function _makeInvisible(node) {
// needs a wrapper around it with contentEditable off to be able to hide it
  node.style.display = "none";
  var el = document.createElement("SPAN");
  el.className = "wrapperForHidden";
  el.contentEditable = false;
  node.parentNode.insertBefore(el,node);
  el.appendChild(node);
}

function _makeVisible(wrapperNode) { 
// needs a wrapper around it with contentEditable off to be able to hide it
  // move the wrapper's children, make them visible, and remove the wrapper
  var p = wrapperNode;
  for (var chld=wrapperNode.firstChild; chld; chld=wrapperNode.firstChild) {
    if (chld.style)
      chld.style.display = "";
    p.parentNode.insertBefore(chld,p);   
  }
  p.parentNode.removeChild(p);
}

HTMLArea.prototype.RpadToggleCodeVisibility = function() {
	var sel = this.getParentElement();
    if (sel.className == "Rpad_input") {
      _makeInvisible(sel);
      return;
    }
    var node = sel;
//alert(node.nodeName);
    while (node.className != "Rpad_results" && node.nodeName != "BODY") {
      node = node.parentNode;
//alert(node.nodeName+' x '+node.innerHTML);
    }
    if (node.nodeName == "BODY") return;
    var previous = node.previousSibling;
    while (previous.nodeType != 1) {
      previous = previous.previousSibling; // skip over text nodes
//alert(previous.nodeName+' x '+previous.innerHTML);
    }
    var isVisible = node.className == "Rpad_results" &&
                    previous.className == "Rpad_input";
    if (isVisible)
      _makeInvisible(previous);
    else if(previous != null && previous.className == "wrapperForHidden")
      _makeVisible(previous);
}

function RpadUnhideElements() {
  var d = document.getElementsByTagName("*");
  for (var i=0; i<=(d.length-1); i++) {
//    d[i].style.display = "";
    if (d[i].className == "wrapperForHidden") 
      _makeVisible(d[i]);
  }
}

function Rpad_calculate_selection() {
  var sel = _Rpad_editor.getParentElement();
  if (sel.className == "Rpad_input") {
    _Rpad_doKeepGoing = false;
    Rpad_calculate_node(sel,true);
  }
  else if (sel.className == "Rpad_results")  {
    _Rpad_doKeepGoing = false;
    Rpad_calculate_node(sel.RpadParent,true);
  }
}

function Rpad_calculate_next(node) {
   _Rpad_doKeepGoing = false;
   Rpad_run_next_node(node, true);
}

function _isInRpadDiv() {
	var parent = _Rpad_editor.getParentElement();

    // find the top level node 2 levels below the BODY
    while (parent.tagName.toLowerCase() != "body") {
      parent = parent.parentNode;
      if (parent.id == 'Rpad_div') return true;
    }
    return false;
}

function Rpad_insert_Rpad(Rpadtype) {
  if (!_isInRpadDiv()) return;
  _Rpad_editor.insertHTML("<pre class='Rpad_input' rpad_type='" + Rpadtype + "'># Enter code or file data here</pre>");
}

function Rpad_insert_Rpad_textarea(Rpadtype) {
  if (!_isInRpadDiv()) return;
  _Rpad_editor.insertHTML("<span contenteditable='false' class='RpadWrapper'><textarea class='Rpad_input' rpad_type='" + Rpadtype + "' rows='5' cols='80'># Enter code or file data here</textarea></span>");
}

function Rpad_insert_Rpad_input(Rpadtype) {
  if (!_isInRpadDiv()) return;
  _Rpad_editor.insertHTML("<span contenteditable='false' class='RpadWrapper'><input class='Rpad_input' rpad_type='" + Rpadtype + "'>&nbsp; </input></span>");
}

function Rpad_insert_Rpad_span(Rpadtype) {
  if (!_isInRpadDiv()) return;
  _Rpad_editor.insertHTML("<span class='Rpad_input' rpad_type='" + Rpadtype + "'>&nbsp; </span>");
}

function Rpad_insert_radio_input() {
  if (!_isInRpadDiv()) return;
  _Rpad_editor.insertHTML('<span contentEditable="false" onMouseDown="this.updateStatusArea"><input type="radio" name="name1" ></input></span>');
}

function Rpad_insert_check_input() {
  if (!_isInRpadDiv()) return;
  _Rpad_editor.insertHTML('<span contentEditable="false"><input type="checkbox" name="name1" ></input></span>');
}
function Rpad_insert_dropdown_input() {
  alert("Not implemented yet");
}

