// Contains replacements for various htmlarea functions.
// Tom Short, tshort@eprisolutions..com 
//
// This is a modified version that replaces the _editorEvent in htmlarea.js.
// The main changes are:
//  * add an F9 key
//  * modify the way "Enter" works on IE
//  * add Ctrl-S for save
//  * take out Ctrl-R (want that for refresh page)
//

// htmlArea v3.0 - Copyright (c) 2002-2004 interactivetools.com, inc.
// This copyright notice MUST stay intact for use (see htmlarea/license.txt).
//
// Portions (c) dynarch.com, 2003-2004
//
// Portions Copyright (c) EPRI Solutions, Inc.., 2004
//
// A free WYSIWYG editor replacement for <textarea> fields.
// For full source code and docs, visit http://www.interactivetools.com/
//
// Version 3.0 developed by Mihai Bazon.
//   http://dynarch.com/mishoo
//

HTMLArea.prototype.isInRpadInput = function() {
	var p = this.getParentElement();
    var result = false;
	while (p && (p.nodeType == 1) && (p.tagName != 'BODY')) {
        result = result || 
                    (p.className == "Rpad_input" && p.tagName != 'INPUT' && p.tagName != 'SELECT' && p.tagName != 'TEXTAREA');
		p = p.parentNode;
	}
	return result;
};

/** A generic event handler for things that happen in the IFRAME's document.
 * This function also handles key bindings. */
HTMLArea.prototype._editorEvent = function(ev) {
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
				range.selectNodeContents(this._doc.body);
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
			this.execCommand(cmd, false, value);
			HTMLArea._stopEvent(ev);
		}
	}
	else if (keyEvent) {
		// other keys here
		switch (ev.keyCode) {
		    case 13: // KEY enter
         		// if (HTMLArea.is_ie) {
                if (this.isInRpadInput()) {
         		  this.insertHTML("<br />");
         		  HTMLArea._stopEvent(ev);
         		}
			break;

		    case 8: // BACKSPACE key - disable for Mozilla/Mozile
                 if (!HTMLArea.is_ie && this.getParentElement().tagName != "INPUT") {
//                   alert(this.getParentElement().tagName);
//                   alert(this.tagName);
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
//                    if (ev.type == "keydown") { // This test distinguishes between F9 & 'x'
                      Rpad_calculate_selection();
				  HTMLArea._stopEvent(ev);
                    }
    			break;
		    case 0x71: // F2 key
                    var isF9 = (!HTMLArea.is_ie || ev.type == "keydown") ;
                    if (isF9) { // This test distinguishes between F9 & 'x'
//                    if (ev.type == "keydown") { // This test distinguishes between F9 & 'x'
                      this.toggleCodeVisibility();
				  HTMLArea._stopEvent(ev);
                    }
    			break;
		}
	}
//editor.updateStatusArea(ev);
	// update the toolbar state after some time
	if (editor._timerToolbar) {
		clearTimeout(editor._timerToolbar);
	}
    var target = HTMLArea.is_ie ? ev.srcElement : ev.target;
	editor._timerToolbar = setTimeout(function() {
		editor.updateStatusArea(target);
		editor.updateToolbar();
		editor._timerToolbar = null;
	}, 50);
};

// retrieve the HTML
HTMLArea.prototype.getHTML = function() {
//    return this._doc.body.innerHTML;
//	return HTMLArea.getHTML(_Rpad_editor._doc.body, true, _Rpad_editor);
//	return HTMLArea.getHTML(this._doc.body, true, this);
	return HTMLArea.getHTML(this._doc.body, false, this);
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

// Replacement for createLink to always make it open in another window
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
		f_target : "_blank"
	};
	this._popupDialog("link.html", function(param) {
		if (!param)
			return false;
		var a = link;
		if (!a) {
			editor._doc.execCommand("createlink", false, param.f_href);
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
		a.target = "_blank";
		a.title = param.f_title.trim();
		editor.selectNodeContents(a);
        // surround the link with a noneditable span so users can activate the link:
        var el = editor._doc.createElement("SPAN");
        el.contentEditable = false;
        el.appendChild(a.cloneNode(true));
        a.replaceNode(el);

		editor.updateToolbar();
	}, outparam);
};
