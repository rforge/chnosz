// Contains various HTMLArea utility functions used in the "lite" editor.
// Tom Short, tshort@eprisolutions..com 
//

// htmlArea v3.0 - Copyright (c) 2002-2004 interactivetools.com, inc.
// This copyright notice MUST stay intact for use (see license.txt).
//
// Portions (c) dynarch.com, 2003-2004
//
// A free WYSIWYG editor replacement for <textarea> fields.
// For full source code and docs, visit http://www.interactivetools.com/
//
// Version 3.0 developed by Mihai Bazon.
//   http://dynarch.com/mishoo
//

var HTMLArea = function() {};

// focuses the iframe window.  returns a reference to the editor document.
HTMLArea.prototype.focusEditor = function() {
  _Rpad.focus(); 
};

/** Returns a node after which we can insert other nodes, in the current
 * selection.  The selection is removed.  It splits a text node, if needed.
 */
HTMLArea.prototype.insertNodeAtSelection = function(toBeInserted) {
	if (!HTMLArea.is_ie) {
		var sel = this._getSelection();
		var range = this._createRange(sel);
		// remove the current selection
		sel.removeAllRanges();
		range.deleteContents();
		var node = range.startContainer;
		var pos = range.startOffset;
		switch (node.nodeType) {
		    case 3: // Node.TEXT_NODE
			// we have to split it at the caret position.
			if (toBeInserted.nodeType == 3) {
				// do optimized insertion
				node.insertData(pos, toBeInserted.data);
				range = this._createRange();
				range.setEnd(node, pos + toBeInserted.length);
				range.setStart(node, pos + toBeInserted.length);
				sel.addRange(range);
			} else {
				node = node.splitText(pos);
				var selnode = toBeInserted;
				if (toBeInserted.nodeType == 11 /* Node.DOCUMENT_FRAGMENT_NODE */) {
					selnode = selnode.firstChild;
				}
				node.parentNode.insertBefore(toBeInserted, node);
				this.selectNodeContents(selnode);
//				this.updateToolbar();
			}
			break;
		    case 1: // Node.ELEMENT_NODE
			var selnode = toBeInserted;
			if (toBeInserted.nodeType == 11 /* Node.DOCUMENT_FRAGMENT_NODE */) {
				selnode = selnode.firstChild;
			}
			node.insertBefore(toBeInserted, node.childNodes[pos]);
			this.selectNodeContents(selnode);
//			this.updateToolbar();
			break;
		}
	} else {
		return null;	// this function not yet used for IE <FIXME>
	}
};

// Returns the deepest node that contains both endpoints of the selection.
HTMLArea.prototype.getParentElement = function() {
	var sel = this._getSelection();
	var range = this._createRange(sel);
	if (HTMLArea.is_ie) {
		switch (sel.type) {
		    case "Text":
		    case "None":
			// It seems that even for selection of type "None",
			// there _is_ a parent element and it's value is not
			// only correct, but very important to us.  MSIE is
			// certainly the buggiest browser in the world and I
			// wonder, God, how can Earth stand it?
			return range.parentElement();
		    case "Control":
			return range.item(0);
		    default:
			return document.body;
		}
	} else try {
		var p = range.commonAncestorContainer;
		if (!range.collapsed && range.startContainer == range.endContainer &&
		    range.startOffset - range.endOffset <= 1 && range.startContainer.hasChildNodes())
			p = range.startContainer.childNodes[range.startOffset];
		/*
		alert(range.startContainer + ":" + range.startOffset + "\n" +
		      range.endContainer + ":" + range.endOffset);
		*/
		while (p.nodeType == 3) {
			p = p.parentNode;
		}
		return p;
	} catch (e) {
		return null;
	}
};

// Returns an array with all the ancestor nodes of the selection.
HTMLArea.prototype.getAllAncestors = function() {
	var p = this.getParentElement();
	var a = [];
	while (p && (p.nodeType == 1) && (p.tagName.toLowerCase() != 'body')) {
		a.push(p);
		p = p.parentNode;
	}
	a.push(document.body);
	return a;
};

// Selects the contents inside the given node
HTMLArea.prototype.selectNodeContents = function(node, pos) {
//	this.forceRedraw();
	var range;
	var collapsed = (typeof pos != "undefined");
	if (HTMLArea.is_ie) {
		range = document.body.createTextRange();
		range.moveToElementText(node);
		(collapsed) && range.collapse(pos);
		range.select();
	} else {
		var sel = this._getSelection();
		range = document.createRange();
		range.selectNodeContents(node);
		(collapsed) && range.collapse(pos);
		sel.removeAllRanges();
		sel.addRange(range);
	}
};

/** Call this function to insert HTML code at the current position.  It deletes
 * the selection, if any.
 */
HTMLArea.prototype.insertHTML = function(html) {
	var sel = this._getSelection();
	var range = this._createRange(sel);
	if (HTMLArea.is_ie) {
		range.pasteHTML(html);
	} else {
		// construct a new document fragment with the given HTML
		var fragment = document.createDocumentFragment();
		var div = document.createElement("div");
		div.innerHTML = html;
		while (div.firstChild) {
			// the following call also removes the node from div
			fragment.appendChild(div.firstChild);
		}
		// this also removes the selection
		var node = this.insertNodeAtSelection(fragment);
	}
};

/**
 *  Call this function to surround the existing HTML code in the selection with
 *  your tags.  FIXME: buggy!  This function will be deprecated "soon".
 */
HTMLArea.prototype.surroundHTML = function(startTag, endTag) {
	var html = this.getSelectedHTML();
	// the following also deletes the selection
	this.insertHTML(startTag + html + endTag);
};

/// Retrieve the selected block
HTMLArea.prototype.getSelectedHTML = function() {
	var sel = this._getSelection();
	var range = this._createRange(sel);
	var existing = null;
	if (HTMLArea.is_ie) {
		existing = range.htmlText;
	} else {
		existing = HTMLArea.getHTML(range.cloneContents(), false, this);
	}
	return existing;
};

/// Return true if we have some selection
HTMLArea.prototype.hasSelectedText = function() {
	// FIXME: come _on_ mishoo, you can do better than this ;-)
	return this.getSelectedHTML() != '';
};

HTMLArea.prototype.stripBaseURL = function(string) {
	var baseurl = document.baseURI || document.URL;
	if (baseurl && baseurl.match(/(.*)\/([^\/]+)/))
		baseurl = RegExp.$1 + "/";

	// strip to last directory in case baseurl points to a file
	baseurl = baseurl.replace(/[^\/]+$/, '');
	var basere = new RegExp(baseurl);
	string = string.replace(basere, "");

	// strip host-part of URL which is added by MSIE to links relative to server root
	baseurl = baseurl.replace(/^(https?:\/\/[^\/]+)(.*)$/, '$1');
	basere = new RegExp(baseurl);
	return string.replace(basere, "");
};

HTMLArea.prototype._createLink = function(link) {
	var editor = this;
	var outparam = null;
	if (typeof link == "undefined") {
		link = this.getParentElement();
		if (link && !/^a$/i.test(link.tagName))
			link = null;
	}
	if (link) outparam = {
		f_href   : HTMLArea.is_ie ? editor.stripBaseURL(link.href) : link.getAttribute("href"),
		f_title  : link.title,
		f_target : link.target
	};
	this._popupDialog("link.html", function(param) {
		if (!param)
			return false;
		var a = link;
		if (!a) {
			document.execCommand("createlink", false, param.f_href);
			a = editor.getParentElement();
			var sel = editor._getSelection();
			var range = editor._createRange(sel);
			if (!HTMLArea.is_ie) {
				a = range.startContainer;
				if (!/^a$/i.test(a.tagName))
					a = a.nextSibling;
			}
		} else a.href = param.f_href.trim();
		if (!/^a$/i.test(a.tagName))
			return false;
		a.target = param.f_target.trim();
		a.title = param.f_title.trim();
		editor.selectNodeContents(a);
//		editor.updateToolbar();
	}, outparam);
};

// retrieve the HTML
HTMLArea.prototype.getHTML = function() {

// need to do this on gecko or textarea changes won't get saved
// note that others are still problematic (inputs)
    if (HTMLArea.is_gecko) {
      var el = _Rpad_editor._doc.getElementsByTagName("textarea");
      for (var i=0; i < el.length; i++)
          el[i].innerHTML = el[i].value;
    }

// Using innerHTML has the following advantages over using DOM traversal:
//      * it's the most complete way to save 
// The advantage of using getHTML (DOM traversal) are:
//      * you can leave out junk that Rpad puts in
//      * the HTML output isn't as junky as innerHTML (on IE anyway)
//      * it doesn't have the IE bug where input elements have "name" attribute fixed

//    return _Rpad.innerHTML;

	return HTMLArea.getHTML(_Rpad, true, this);
};

// retrieve the HTML (fastest version, but uses innerHTML)
HTMLArea.prototype.getInnerHTML = function() {
	switch (this._editMode) {
	    case "wysiwyg"  :
		if (!this.config.fullPage)
			return this._doc.body.innerHTML;
		else
			return this.doctype + "\n" + this._doc.documentElement.innerHTML;
	    case "textmode" : return this._textArea.value;
	    default	    : alert("Mode <" + mode + "> not defined!");
	}
	return false;
};

// completely change the HTML inside
HTMLArea.prototype.setHTML = function(html) {
	switch (this._editMode) {
	    case "wysiwyg"  :
		if (!this.config.fullPage)
			this._doc.body.innerHTML = html;
		else
			// this._doc.documentElement.innerHTML = html;
			this._doc.body.innerHTML = html;
		break;
	    case "textmode" : this._textArea.value = html; break;
	    default	    : alert("Mode <" + mode + "> not defined!");
	}
	return false;
};

// sets the given doctype (useful when config.fullPage is true)
HTMLArea.prototype.setDoctype = function(doctype) {
	this.doctype = doctype;
};



/***************************************************
 *  Category: UTILITY FUNCTIONS
 ***************************************************/

// browser identification

HTMLArea.agt = navigator.userAgent.toLowerCase();
HTMLArea.is_ie	   = ((HTMLArea.agt.indexOf("msie") != -1) && (HTMLArea.agt.indexOf("opera") == -1));
HTMLArea.is_opera  = (HTMLArea.agt.indexOf("opera") != -1);
HTMLArea.is_mac	   = (HTMLArea.agt.indexOf("mac") != -1);
HTMLArea.is_mac_ie = (HTMLArea.is_ie && HTMLArea.is_mac);
HTMLArea.is_win_ie = (HTMLArea.is_ie && !HTMLArea.is_mac);
HTMLArea.is_gecko  = (navigator.product == "Gecko");

// variable used to pass the object to the popup editor window.
HTMLArea._object = null;

// function that returns a clone of the given object
HTMLArea.cloneObject = function(obj) {
	var newObj = new Object;

	// check for array objects
	if (obj.constructor.toString().indexOf("function Array(") == 1) {
		newObj = obj.constructor();
	}

	// check for function objects (as usual, IE is fucked up)
	if (obj.constructor.toString().indexOf("function Function(") == 1) {
		newObj = obj; // just copy reference to it
	} else for (var n in obj) {
		var node = obj[n];
		if (typeof node == 'object') { newObj[n] = HTMLArea.cloneObject(node); }
		else                         { newObj[n] = node; }
	}

	return newObj;
};

// FIXME!!! this should return false for IE < 5.5
HTMLArea.checkSupportedBrowser = function() {
	if (HTMLArea.is_gecko) {
		if (navigator.productSub < 20021201) {
			alert("You need at least Mozilla-1.3 Alpha.\n" +
			      "Sorry, your Gecko is not supported.");
			return false;
		}
		if (navigator.productSub < 20030210) {
			alert("Mozilla < 1.3 Beta is not supported!\n" +
			      "I'll try, though, but it might not work.");
		}
	}
	return HTMLArea.is_gecko || HTMLArea.is_ie;
};

// selection & ranges

// returns the current selection object
HTMLArea.prototype._getSelection = function() {
	if (HTMLArea.is_ie) {
		return document.selection;
	} else {
		return window.getSelection();
	}
};

// returns a range for the current selection
HTMLArea.prototype._createRange = function(sel) {
	if (HTMLArea.is_ie) {
		return sel.createRange();
	} else {
		if (typeof sel != "undefined") {
			try {
				return sel.getRangeAt(0);
			} catch(e) {
				return document.createRange();
			}
		} else {
			return document.createRange();
		}
	}
};

// event handling

HTMLArea._addEvent = function(el, evname, func) {
	if (HTMLArea.is_ie) {
		el.attachEvent("on" + evname, func);
	} else {
		el.addEventListener(evname, func, true);
	}
};

HTMLArea._addEvents = function(el, evs, func) {
	for (var i in evs) {
		HTMLArea._addEvent(el, evs[i], func);
	}
};

HTMLArea._removeEvent = function(el, evname, func) {
	if (HTMLArea.is_ie) {
		el.detachEvent("on" + evname, func);
	} else {
		el.removeEventListener(evname, func, true);
	}
};

HTMLArea._removeEvents = function(el, evs, func) {
	for (var i in evs) {
		HTMLArea._removeEvent(el, evs[i], func);
	}
};

HTMLArea._stopEvent = function(ev) {
	if (HTMLArea.is_ie) {
		ev.cancelBubble = true;
		ev.returnValue = false;
	} else {
		ev.preventDefault();
		ev.stopPropagation();
	}
};

HTMLArea._removeClass = function(el, className) {
	if (!(el && el.className)) {
		return;
	}
	var cls = el.className.split(" ");
	var ar = new Array();
	for (var i = cls.length; i > 0;) {
		if (cls[--i] != className) {
			ar[ar.length] = cls[i];
		}
	}
	el.className = ar.join(" ");
};

HTMLArea._addClass = function(el, className) {
	// remove the class first, if already there
	HTMLArea._removeClass(el, className);
	el.className += " " + className;
};

HTMLArea._hasClass = function(el, className) {
	if (!(el && el.className)) {
		return false;
	}
	var cls = el.className.split(" ");
	for (var i = cls.length; i > 0;) {
		if (cls[--i] == className) {
			return true;
		}
	}
	return false;
};

HTMLArea.isBlockElement = function(el) {
	var blockTags = " body form textarea fieldset ul ol dl li div " +
		"p h1 h2 h3 h4 h5 h6 quote pre table thead " +
		"tbody tfoot tr td iframe address ";
	return (blockTags.indexOf(" " + el.tagName.toLowerCase() + " ") != -1);
};

HTMLArea.needsClosingTag = function(el) {
	var closingTags = " head script style div span tr td tbody table em strong font a title ";
	return (closingTags.indexOf(" " + el.tagName.toLowerCase() + " ") != -1);
};

// performs HTML encoding of some given string
HTMLArea.htmlEncode = function(str) {
	// we don't need regexp for that, but.. so be it for now.
	str = str.replace(/&/ig, "&amp;");
	str = str.replace(/</ig, "&lt;");
	str = str.replace(/>/ig, "&gt;");
	str = str.replace(/\x22/ig, "&quot;");
	// \x22 means '"' -- we use hex reprezentation so that we don't disturb
	// JS compressors (well, at least mine fails.. ;)
	return str;
};

// Replaces the version in htmlarea.js to not "avoid certain attributes".
//
HTMLArea.getHTML = function(root, outputRoot, editor) {
	var html = "";
	switch (root.nodeType) {
	    case 1: // Node.ELEMENT_NODE
	    case 11: // Node.DOCUMENT_FRAGMENT_NODE
		var closed;
		var i;
		var root_tag = (root.nodeType == 1) ? root.tagName.toLowerCase() : '';
		if (HTMLArea.is_ie && root_tag == "head") {
			if (outputRoot)
				html += "<head>";
			// lowercasize
			var save_multiline = RegExp.multiline;
			RegExp.multiline = true;
			var txt = root.innerHTML.replace(HTMLArea.RE_tagName, function(str, p1, p2) {
				return p1 + p2.toLowerCase();
			});
			RegExp.multiline = save_multiline;
			html += txt;
			if (outputRoot)
				html += "</head>";
			break;
		} else if (outputRoot) {
			closed = (!(root.hasChildNodes() || HTMLArea.needsClosingTag(root)));
			html = "<" + root.tagName.toLowerCase();
			var attrs = root.attributes;
			for (i = 0; i < attrs.length; ++i) {
				var a = attrs.item(i);
//				if (!a.specified) {
//bugfix: http://www.interactivetools.com/forum/forum.cgi?do=post_view_printable;post=23072;guest=5256845&t=search_engine
                if (!a.specified && !(a.nodeName.toLowerCase() == "value" && a.nodeValue != "")
                                 && !(a.nodeName.toLowerCase() == "name" && a.nodeValue != "")) {
					continue;
				}
				var name = a.nodeName.toLowerCase();
				if (/_moz|contenteditable/.test(name)) {
					// avoid certain attributes
// Don't avoid these attributes!!
//					continue;
				}
				var value;
				if (name != "style") {
					// IE5.5 reports 25 when cellSpacing is
					// 1; other values might be doomed too.
					// For this reason we extract the
					// values directly from the root node.
					// I'm starting to HATE JavaScript
					// development.  Browser differences
					// suck.
					//
					// Using Gecko the values of href and src are converted to absolute links
					// unless we get them using nodeValue()
					if (typeof root[a.nodeName] != "undefined" && name != "href" && name != "src") {
//						value = root[a.nodeName];
						value = a.nodeValue;
					} else {
						value = a.nodeValue;
						// IE seems not willing to return the original values - it converts to absolute
						// links using a.nodeValue, a.value, a.stringValue, root.getAttribute("href")
						// So we have to strip the baseurl manually -/
						if (HTMLArea.is_ie && (name == "href" || name == "src")) {
							value = editor.stripBaseURL(value);
						}
					}
				} else { // IE fails to put style in attributes list
					// FIXME: cssText reported by IE is UPPERCASE
					value = root.style.cssText;
				}
				if (/(_moz|^$)/.test(value)) {
					// Mozilla reports some special tags
					// here; we don't need them.
					continue;
				}
				html += " " + name + '="' + value + '"';
			}
			html += closed ? " />" : ">";
		}
		if (HTMLArea.is_ie && outputRoot && !closed && root.tagName.toLowerCase() == "script") {
            html += root.text;
        }
		for (i = root.firstChild; i; i = i.nextSibling) {
			html += HTMLArea.getHTML(i, true, editor);
		}
		if (outputRoot && !closed) {
			html += "</" + root.tagName.toLowerCase() + ">";
		}
		break;
	    case 3: // Node.TEXT_NODE
		html = HTMLArea.htmlEncode(root.data);
		break;
	    case 8: // Node.COMMENT_NODE
		html = "<!--" + root.data + "-->";
		break;		// skip comments, for now.
	}
	return html;
};

/**
 * FIX: Internet Explorer returns an item having the _name_ equal to the given
 * id, even if it's not having any id.  This way it can return a different form
 * field even if it's not a textarea.  This workarounds the problem by
 * specifically looking to search only elements having a certain tag name.
 */
HTMLArea.getElementById = function(tag, id) {
	var el, i, objs = document.getElementsByTagName(tag);
	for (i = objs.length; --i >= 0 && (el = objs[i]);)
		if (el.id == id)
			return el;
	return null;
};



// EOF
// Local variables: //
// c-basic-offset:8 //
// indent-tabs-mode:t //
// End: //
