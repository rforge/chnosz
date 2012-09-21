/*
    Rpad_editor.js  --  creates the "full" editor based on HTMLArea that works in IE or Mozilla

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

// Add a textarea wrapper around the whole document body.
// This is somewhat cludgy, but htmlarea is a replacement for textarea's and not a general page editor.
// Also, watch out that the wrapper must skip over the first document element, which should be the table
// that holds the menu bar.
function add_textarea_wrapper() {
    var d = document.getElementsByTagName("body")[0];
    var textarea = document.createElement('textarea');
    textarea.id = "Rpad";
    textarea.value = d.innerHTML;
    d.insertBefore(textarea,d.firstChild);
    while (d.firstChild.nextSibling) // zap the body's children
      d.removeChild(d.firstChild.nextSibling);
}  
   
// This is a modified version of the resize_editor in htmlarea.js used for fullscreen operation.
function resize_editor() {  // resize editor to fix window
  var newHeight;
  if (document.all) {
    // IE
    newHeight = document.body.offsetHeight - _Rpad_editor._toolbar.offsetHeight 
                - _Rpad_editor._menuBar.offsetHeight - 25;
  } else {
    // Gecko
    newHeight = window.innerHeight - _Rpad_editor._toolbar.offsetHeight 
                - _Rpad_editor._menuBar.offsetHeight - 20;
  }
  if (_Rpad_editor.config.statusBar) {
    newHeight -= _Rpad_editor._statusBar.offsetHeight;
  }
  if (newHeight < 0) { newHeight = 50; }
  _Rpad_editor._textArea.style.height = _Rpad_editor._iframe.style.height = newHeight+'px';
}  
   
var _Rpad_editor = null;
var _Rpad = null;

// This is the main initialization routine.
function initEditor() {
  add_textarea_wrapper();
  _Rpad_editor = new HTMLArea("Rpad");
  var cfg = _Rpad_editor.config; // this is the default configuration
  cfg.width	   = "100%";
  cfg.height	   = "auto";
  cfg.toolbar = [                                                                                                                                                                        
    	  [                                                                                                                                                                              
    	  "bold", "italic", "underline", "separator",                                                                                                                                    
    	  "strikethrough", "subscript", "superscript", "separator",                                                                                                                      
//  	  "copy", "cut", "paste", "space",                                                                                                                                               
          "undo", "redo",                                                                                                                                                                
          "justifyleft", "justifycenter", "justifyright", "justifyfull", "separator",                                                                                                    
    	  "insertorderedlist", "insertunorderedlist", "outdent", "indent", "separator",                                                                                                  
    	  "forecolor", "hilitecolor", "separator",                                                                                                                                       
    	  "inserthorizontalrule", "createlink", "insertimage", "inserttable", "htmlmode",                                                                                                
          "fontname", "space",                                                                                                                                                           
    	  "fontsize", "space",                                                                                                                                                           
    	  "formatblock", "space"                                                                                                                                                         
          ]                                                                                                                                                                              
    	];                                                                                                                                                                               
  _Rpad_editor._customUndo = true;
                                                                                                                                                                                         
// hardcode this stuff in because I can't get Mozilla to read the imported css file using @import                                                                                                     
  cfg.pageStyle = "body { background-color: #e4eeee; font-family: tahoma; font-size: 18px; padding:5; margin:5}"+                                                                        
                  ".Rpad_input "+                                                                                                                                                        
                       "{ background-color: white; font-family: monospace;} "+                                                                                                           
                  ".Rpad_results "+                                                                                                                                                      
                       "{ background-color: infobackground; font-family: monospace;} "+                                                                                                  
"H1 { font-family: Arial, Helvetica, sans-serif; font-size: 20pt; font-style: normal; font-weight: normal; color: #FFFFFF; background: #004080; text-align: center; margin: 10pt 2.5%}"+  
"H2 { font-family: Arial, Helvetica, sans-serif; font-size: 18pt; font-style: normal; font-weight: normal; color: #FFFFFF; background: #000000}                                       "+  
"H2.index { font-family: Arial, Helvetica, sans-serif; font-size: 18pt; font-style: normal; font-weight: normal; color: #FFFFFF; background: #000000; margin: 10pt 5%}                "+ 
"H3 { font-family: Arial, Helvetica, sans-serif; font-size: 14pt; font-style: normal; font-weight: bold; color: #004080  }                                                            "+ 
"H4 { font-family: T, Helvetica, sans-serif; font-size: 10pt; font-style: normal; font-weight: bold; color: #000000; line-height: 16pt}                                               "+ 
"LI { font-family: Verdana, Arial, Helvetica, sans-serif; font-size: 10pt }                                                                                                           "+ 
"A { font-family: Verdana, Arial, Helvetica, sans-serif; font-size: 10pt ;text-decoration:none}                                                                                         "+
".caption{font-style:italic}                                                                                                                                                            "+
".title2{ font-family: Arial, Helvetica, sans-serif; font-size: 14pt; font-style: normal; font-weight: bold; color: #004080  }                                                          "+
"P { font-family: Verdana, Arial, Helvetica, sans-serif; font-size: 10pt ;font-weight:normal}                                                                                           "+
".command {font-family=verdana, arial; color=red	}                                                                                                                                   "+
".function{font-family=verdana, arial; color=blue}                                                                                                                                      "+
".partitle{font-family=verdana, arial; font-weight:bold}                                                                                                                                "+
"XMP{font-family: Verdana, Arial, Helvetica, sans-serif; font-size: 10pt }                                                                                                              "+
".function{font-family=courier; color=blue; font-size: 10pt }                                                                                                                           "+
"TABLE.dataframe{align:left; background: #FFFFFF}                                                                                                                                       "+
"TR{font-family: Arial, Helvetica,Times, Helvetica, sans-serif; font-size: 10pt; font-style: normal  ; padding:0 0}                                                                     "+
"TR.firstline{color: #FFFFFF; background: #000000;text-align=center;font-weight: bold}                                                                                                  "+
"TR.ListBackTitle{color: #FFFFFF; background: #000000;text-align=left;font-weight: bold}                                                                                                "+
"TD{padding:0 0}                                                                                                                                                                        "+
"TD.ListBackMain{background: #E0E0E0; padding:0 0}                                                                                                                                      "+
"TD.firstcolumn{padding:5 10;background:#C0C0C0;text-align=right}                                                                                                                       "+
"TD.cellinside{padding:5 10;background:#FFFFFF;text-align=right}                                                                                                                        "+
                  "";                                                                                                                                                                  
  cfg.pageStyle += "\n@import url(editor/Rpad_page.css);"; // this doesn't seem to do anything under Mozilla
                                                                                                                                                                                       
  _Rpad_editor.generate();                                                                                                                                                             
                                                                                                                                                                                       
                                                                                                                                                                                       
                                                                                                                                                                                       
  // set up the menu bar and status area                                                                                                                                               
                                                                                                                                                                                       
  var el = document.createElement("div");                                                                                                                                              
  el.innerHTML = '<table height="30px" width="100%"><td><div id="Rpad_menuBar" unselectable="on"></div></td><td align="right"><div id="Rpad_statusArea" style="fontsize:3"></div></td></table>';
  document.body.insertBefore(el,document.body.firstChild);                                                                                                                             
                                                                                                                                                                                       
  cmDraw('Rpad_menuBar', Rpad_menuBar, 'hbr', cmThemeOffice, 'ThemeOffice');                                                                                                           
  _Rpad_editor._menuBar = document.getElementById("Rpad_menuBar");                                                                                                                     
  _Rpad_editor._statusArea = document.getElementById("Rpad_statusArea");
  _Rpad_editor._iframe.style.width = "100%";
  _Rpad_editor._textArea.style.width = "100%";


  resize_editor();

  // set child window contents and event handlers, after a small delay
  setTimeout(function() {
    // setup event handlers
    _Rpad = _Rpad_editor._doc.getElementsByTagName("body")[0];
   
    window.onresize = resize_editor;
    // HTMLArea modifies onunload internally, so we want to change it back.
    window.onunload = stopRpad;
  }, 333);                      // give it some time to meet the new frame
}

HTMLArea.prototype._createMenuBar = function() {
    var el = document.createElement("span");
    el.id = "Rpad_menuBar1";
    el.className = "menuBar";  
    el.style.fontSize = "12";

    this._htmlArea.appendChild(el);
    this._menuBar = el;
};


HTMLArea.prototype.updateStatusArea = function(node) {
  var parent = this.getParentElement();
  if (node.tagName == 'SELECT' || node.tagName == 'INPUT') // try an INPUT or SELECT
    parent = node;

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
    inp.value = parent.getAttribute("name");
    inp.onchange = function() {parent.setAttribute("name",this.value);parent.name = this.value;};
    inp.onblur = function() {parent.setAttribute("name",this.value);parent.name = this.value;};
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
   		[null, 'Unhide all code', "javascript:_Rpad_editor.unhideElements()", '_self', null],
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
  var el = _Rpad_editor._doc.createElement("SPAN");
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

HTMLArea.prototype.toggleCodeVisibility = function() {
    this.focusEditor();
	var sel = this.getParentElement();
    if (sel.className == "Rpad_input") {
      _makeInvisible(sel);
      return;
    }
    var previous = sel.previousSibling;
    while (previous.nodeType != 1) previous = previous.previousSibling; // skip over text nodes

    var isVisible = sel.className == "Rpad_results" &&
                    previous.className == "Rpad_input";
    if (isVisible)
      _makeInvisible(previous);
    else if(previous != null && previous.className == "wrapperForHidden")
      _makeVisible(previous);
}

HTMLArea.prototype.unhideElements = function() {
  this.focusEditor();
  var d = _Rpad_editor._doc.getElementsByTagName("*");
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

function Rpad_insert_Rpad(Rpadtype) {
  _Rpad_editor.focusEditor();
  _Rpad_editor.insertHTML("<pre class='Rpad_input' rpad_type='" + Rpadtype + "'># Enter code or file data here</pre>");
}

function Rpad_insert_Rpad_textarea(Rpadtype) {
  _Rpad_editor.focusEditor();
  _Rpad_editor.insertHTML("<span contenteditable='false' class='RpadWrapper'><textarea class='Rpad_input' rpad_type='" + Rpadtype + "' rows='5' cols='80'># Enter code or file data here</textarea></span>");
}

function Rpad_insert_Rpad_input(Rpadtype) {
  _Rpad_editor.focusEditor();
  _Rpad_editor.insertHTML("<span contenteditable='false'><input class='Rpad_input' rpad_type='" + Rpadtype + "'>&nbsp; </input></span>");
}

function Rpad_insert_Rpad_span(Rpadtype) {
  _Rpad_editor.focusEditor();
  _Rpad_editor.insertHTML("<span class='Rpad_input' rpad_type='" + Rpadtype + "'>&nbsp; </span>");
}

function Rpad_insert_radio_input() {
  _Rpad_editor.focusEditor();
  _Rpad_editor.insertHTML('<span contentEditable="false" onMouseDown="this.updateStatusArea"><input type="radio" name="name1" ></input></span>');
}

function Rpad_insert_check_input() {
  _Rpad_editor.focusEditor();
  _Rpad_editor.insertHTML('<span contentEditable="false"><input type="checkbox" name="name1" ></input></span>');
}
function Rpad_insert_dropdown_input() {
  alert("Not implemented yet");
}

