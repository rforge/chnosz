/*
    Rpad_loader.js  --  loads the necessary javascript files
                    --  change "Rpad_editmode" to switch between the "lite" and "heavy" editor

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

if (Rpad_editmode == null)  var Rpad_editmode = "default";  // "lite" for lite editor, "heavy" for HTMLArea editor, other for defaults
//var Rpad_editmode = "lite";
var is_gecko  = (navigator.product == "Gecko");

// default is "lite" for gecko (mozilla) and "heavy" for IE

if (Rpad_editmode == "lite" || (is_gecko && Rpad_editmode != "heavy")) {  

  document.write('<script src="editor/JSCookMenu/JSCookMenu.js" rpadIgnore="true"></script>');
  document.write('<link rel="stylesheet" href="editor/JSCookMenu/ThemeOffice/theme.css" type="text/css" rpadIgnore="true">');
  document.write('<script src="editor/JSCookMenu/ThemeOffice/theme.js" rpadIgnore="true"></script>');
   
  document.write('<style type="text/css" rpadIgnore="true">@import url(editor/Rpad.css);</style>');
  document.write('<script src="editor/htmlarea_utils.js" rpadIgnore="true"> </script>');
  document.write('<script src="editor/Rpad_processing_asynch.js" rpadIgnore="true"  > </script>');
  document.write('<script src="editor/Rpad_lite_editor.js" rpadIgnore="true" > </script>');
 
} else { // heavy

  document.write('<script src="editor/JSCookMenu/JSCookMenu.js" rpadIgnore="true"></script>');
  document.write('<link rel="stylesheet" href="editor/JSCookMenu/ThemeOffice/theme.css" type="text/css" rpadIgnore="true">');
  document.write('<script src="editor/JSCookMenu/ThemeOffice/theme.js" rpadIgnore="true"></script>');
      
  document.write('<style type="text/css" rpadIgnore="true">@import url(editor/htmlarea/htmlarea.css);</style>');
  document.write('<style type="text/css" rpadIgnore="true">@import url(editor/Rpad_wrapper.css);</style>');
  document.write('<script type="text/javascript" rpadIgnore="true">_editor_url = "editor/htmlarea"; _editor_lang = "en";</script>');
  document.write('<script src="editor/htmlarea/htmlarea.js" rpadIgnore="true"> </script>');
  document.write('<script src="editor/htmlarea/lang/en.js" rpadIgnore="true" > </script>');
  document.write('<script src="editor/htmlarea/dialog.js" rpadIgnore="true"  > </script>');
  
  document.write('<script src="editor/Rpad_processing_asynch.js" rpadIgnore="true"  > </script>');
  document.write('<script src="editor/htmlarea_replacements.js" rpadIgnore="true" > </script>');
  document.write('<script src="editor/Rpad_editor.js" rpadIgnore="true" > </script>');

}