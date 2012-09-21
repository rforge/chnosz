/*
    Rpad_processing.js  --  handles communication with the server (asynchronously), 
                            parsing of the HTML page, and displaying of results

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


function myXmlHttp() {
// from: http://jibbering.com/2002/4/httprequest.html
  var xmlhttp
  /*@cc_on @*/
  /*@if (@_jscript_version >= 5)
    try {
    xmlhttp=new ActiveXObject("Msxml2.XMLHTTP")
   } catch (e) {
    try {
      xmlhttp=new ActiveXObject("Microsoft.XMLHTTP")
    } catch (E) {
     xmlhttp=false
    }
   }
  @else
   xmlhttp=false
  @end @*/
  if (!xmlhttp) {
   try {
//    netscape.security.PrivilegeManager.enablePrivilege("UniversalBrowserRead");
    xmlhttp = new XMLHttpRequest();
   } catch (e) {
    xmlhttp=false
   }
  }
  return xmlhttp;
}
// end from: http://jibbering.com

var Rpad_dir = "";
var _Rpad_Results = "";
var _runState = "init";  

var _RProcessStr = "server/R_process.pl";
if (document.URL.search(/^file:/) >= 0) { // the file is being run locally (no server)
  _RProcessStr = "http://localhost:80";
  var _isLocal = true;
  var _isLocalAndMoz = HTMLArea.is_gecko; // we already know it's local
  if (_isLocalAndMoz) netscape.security.PrivilegeManager.enablePrivilege('UniversalBrowserRead');
} 
// setwd('e:/misc/www/Rpadl/server/local')
// source('R_server.R')
// startRpadServer()



function Rpad_login() {
  if (_isLocal) { R_login(); Rpad_dir='local'; return;}
  var xmlhttp = myXmlHttp();
  xmlhttp.open("POST", "server/Rpad_process.pl",true); //asynchronous connection
  xmlhttp.setRequestHeader("Content-Type", "application/x-www-form-urlencoded"); 
  xmlhttp.onreadystatechange=function() {
    if (xmlhttp.readyState==4) {
      Rpad_dir = xmlhttp.responseText;
      window.status = 'Rpad logged in';
      R_login();
    }
  }
  xmlhttp.send("command=login"); 
  window.status = 'Rpad logging in';
}

function Rpad_logout() {
  if (_isLocal) return;
  var xmlhttp = myXmlHttp();
  if (_isLocalAndMoz) netscape.security.PrivilegeManager.enablePrivilege("UniversalBrowserRead");
  xmlhttp.open("POST", "server/Rpad_process.pl",true); //asynchronous connection - send & go
  xmlhttp.setRequestHeader("Content-Type", "application/x-www-form-urlencoded"); 
  xmlhttp.send("command=logout&ID="+Rpad_dir); 
}

function R_login() {
  if (_isLocal) return;
  var xmlhttp = myXmlHttp();
  if (_isLocalAndMoz) netscape.security.PrivilegeManager.enablePrivilege("UniversalBrowserRead");
  xmlhttp.open("POST", _RProcessStr, true); //asynchronous connection
  xmlhttp.setRequestHeader("Content-Type", "application/x-www-form-urlencoded"); 
  xmlhttp.send("ID="+Rpad_dir+"&command=login"); 
}

function R_logout() {
  if (_isLocal) return;
  var xmlhttp = myXmlHttp();
  if (_isLocalAndMoz) netscape.security.PrivilegeManager.enablePrivilege("UniversalBrowserRead");
  xmlhttp.open("POST", _RProcessStr,true); //asynchronous connection
  xmlhttp.setRequestHeader("Content-Type", "application/x-www-form-urlencoded"); 
  xmlhttp.onreadystatechange=function() {
    if (xmlhttp.readyState==4) {
      Rpad_logout();
    }
  }
  xmlhttp.send("ID="+Rpad_dir+"&command=logout"); 
//with the following alert, Mozilla will terminate the R process
//without the alert, it gives an error on the server and the R process stays alive
//alert('here3b: '+Rpad_dir);
}

function fixhtmlencodings(str) {  
  str = str.replace(/&/g,"&amp;");  
  str = str.replace(/</g,"&lt;");  
  str = str.replace(/>/g,"&gt;");  
  str = str.replace(/\n/g,"\n<BR/>");
  str = str.replace(/ /g,"&#160;"); // replace spaces with nonbreaking spaces
  return(str);
}

function fixformat(str,type) {
// <htmlon/> turns html mode on
// <htmloff/> turns html mode off
// s='aasdf asdf asdf <htmlon/>asdf asdf <htmloff/>asdf asdf'
  var resultstr = "";
  var ishtml = type == 'html';
  var s = "";

  do {
    if (ishtml) {
      var s1 = str.split('<htmloff\/>');
      s1[0] = s1[0].replace(/<htmlon\/>/gi,""); // get rid of redundant tags
      s = s1[0];
      s1[0] = "";
      str = s1.join('<htmloff\/>').replace('<htmloff\/>','') // rejoin all but the first with <htmloff/>
      resultstr += s;
      ishtml = false;
    } else {
      var s1 = str.split('<htmlon\/>');
      s1[0] = s1[0].replace(/<htmloff\/>/gi,""); // get rid of redundant tags
      s = s1[0];
      s1[0] = "";
      str = s1.join('<htmlon\/>').replace('<htmlon\/>','')
      resultstr += fixhtmlencodings(s);
      ishtml = true;
    }
  } while (str != "");
  return(resultstr);
}

    
function trim(str) {
  return(str.replace(/^(\s)(\w\W)(\s)$/g,"$2"));  
}

var _Rpad_doKeepGoing = true;
var _Rpad_hasResults = true;
var _serverTries = 0;  

function R_run_node(commands,node) {
  var xmlhttp = myXmlHttp();
  if (_isLocalAndMoz) netscape.security.PrivilegeManager.enablePrivilege("UniversalBrowserRead");
  window.status = 'Sending commands to the R calculation server';
  if(Rpad_dir == "") {
    _serverTries = _serverTries + 1;
    window.status = 'Sending commands to the R calculation server: try ' + _serverTries;
    if (_serverTries > 5)
    alert("Can't connect to the server");
    else
      setTimeout(function() {
        R_run_node(commands,node);
      }, 2500);  // give it some time & try again
    return;
  }
  xmlhttp.open("POST", _RProcessStr, true);
  xmlhttp.setRequestHeader("Content-Type", "application/x-www-form-urlencoded"); 
  xmlhttp.onreadystatechange=function() {
    if (xmlhttp.readyState==4) {
      var result = fixformat(trim(xmlhttp.responseText),_Rpad_Results);
      window.status='R results received';
      if (_Rpad_Results == 'none') result = "";
      append_results(result,node);
      if (_Rpad_doKeepGoing)
        Rpad_run_next_node(node);
    }
  }
  xmlhttp.send("ID="+Rpad_dir+"&command=R_commands&R_commands="+encodeURIComponent(commands)); 
}

function R_run_commands(commands, jsOnDone) {
  var xmlhttp = myXmlHttp();
  if (_isLocalAndMoz) netscape.security.PrivilegeManager.enablePrivilege("UniversalBrowserRead");
  window.status = 'Sending commands to the R calculation server';
  if(Rpad_dir == "") {
    _serverTries = _serverTries + 1;
    window.status = 'Sending commands to the R calculation server: try ' + _serverTries;
    if (_serverTries > 5)
    alert("Can't connect to the server");
    else
      setTimeout(function() {
        R_run_node(commands,node);
      }, 2500);  // give it some time & try again
    return;
  }
  xmlhttp.open("POST", _RProcessStr, true);
  xmlhttp.onreadystatechange=function() {
    if (xmlhttp.readyState==4) {
      var result = fixformat(trim(xmlhttp.responseText),_Rpad_Results);
      window.status='R results received';
      if (jsOnDone != null)
        eval(jsOnDone);
    }
  }
  xmlhttp.setRequestHeader("Content-Type", "application/x-www-form-urlencoded"); 
  xmlhttp.send("ID="+Rpad_dir+"&command=R_commands&R_commands="+encodeURIComponent(commands)); 
}

function shell_run_node(commands,node) {
  var xmlhttp = myXmlHttp();
  if(Rpad_dir == "") {
    alert("Can't connect to the server");
    return;
  }
  window.status = 'Sending shell commands to the server';
  xmlhttp.open("POST", "server/shell_process.pl",true);
  xmlhttp.setRequestHeader("Content-Type", "application/x-www-form-urlencoded"); 
  xmlhttp.onreadystatechange=function() {
    if (xmlhttp.readyState==4) {
      var result = fixformat(trim(xmlhttp.responseText),_Rpad_Results);
      window.status = 'shell results received';
      if (_Rpad_Results == 'none') result = "";
      append_results(result,node);
      if (_Rpad_doKeepGoing)
        Rpad_run_next_node(node);
    }
  }
  xmlhttp.send("ID="+Rpad_dir+"&shell_commands="+encodeURIComponent(commands)); 
}

function put_file_node(contents,filename,node) {
  var xmlhttp = myXmlHttp();
  if(Rpad_dir == "") {
    alert("Can't connect to the server");
    return;
  }
  window.status = 'Sending file to the server';
  xmlhttp.open("POST", "server/Rpad_process.pl",true);
  xmlhttp.onreadystatechange=function() {
    if (xmlhttp.readyState==4) {
      var result = fixformat(trim(xmlhttp.responseText));
      window.status = 'File sent';
      if (_Rpad_doKeepGoing)
        Rpad_run_next_node(node);
    }
  }
  xmlhttp.setRequestHeader("Content-Type", "application/x-www-form-urlencoded"); 
  var s = "ID="+Rpad_dir+"&command=savefile&filename=" + filename + "&content=" + encodeURIComponent(contents); 
  xmlhttp.send(s); 
}

function append_results(result,node) {
  if (result=="") {
    var hasResults = node.nextSibling != null && node.nextSibling.className == "Rpad_results";
    if (hasResults)
      node.parentNode.removeChild(node.nextSibling);
    return;
  }
  var isBlock = node.nodeName == "DIV" || node.nodeName == "P" || node.nodeName == "TEXTAREA" ||
                node.nodeName == "PRE";
  if (isBlock) {
    var resultNode = _Rpad_editor._doc.createElement("DIV");
  } else {
    var resultNode = _Rpad_editor._doc.createElement("SPAN");
  }
  resultNode.className = "Rpad_results";
  resultNode.style.display = "";
  resultNode.innerHTML = result;
  resultNode.RpadParent = node;
  var parentHasWrapper = node.parentNode.className == "wrapperForHidden" ||
                         node.parentNode.className == "RpadWrapper";
  if (parentHasWrapper)
    node = node.parentNode;
  var missingResults = node.nextSibling == null || node.nextSibling.className != "Rpad_results";
  if (missingResults)
    node.parentNode.insertBefore(resultNode,node.nextSibling);
  else
    node.parentNode.replaceChild(resultNode,node.nextSibling);
}

function append_results_node(newnode,node) {
  var resultNode = _Rpad_editor._doc.createElement(node.tagName);
  resultNode.className = "Rpad_results";
  resultNode.style.display = "";
  resultNode.appendChild(newnode);
  resultNode.RpadParent = node;
  var parentHasWrapper = node.parentNode.className == "wrapperForHidden";
  if (parentHasWrapper)
    node = node.parentNode;
  var missingResults = node.nextSibling == null || node.nextSibling.className != "Rpad_results";
  if (missingResults)
    node.parentNode.insertBefore(resultNode,node.nextSibling);
  else
    node.parentNode.replaceChild(resultNode,node.nextSibling);
}
   
function getTextContent(node) {
  if (node.nodeType == 3) { // text node
    if (node.parentNode.nodeName == "PRE") { // really should see if the style has "white-space: pre"
      return node.nodeValue;
    } else {
      return node.nodeValue.replace(/\n/g," ");
    }
  }
  if (node.nodeType == 1) { // element node
    if (node.nodeName == "BR") return "\n";
//    if (node.nodeName == "TEXTAREA") return node.innerHTML;
    if (node.nodeName == "TEXTAREA") return node.value;
    if (node.nodeName == "INPUT") return node.value;
    var text = [];
    for (var chld = node.firstChild;chld;chld=chld.nextSibling) {
        text.push(getTextContent(chld));
    } return text.join("");
  } return ""; // some other node, won't contain text nodes.
}  

function Rpad_calculate() {
//  var node = _Rpad_editor._doc.getElementsByTagName("body")[0];
  var node = _Rpad;
  _Rpad_doKeepGoing = true;
  _serverTries = 0;
  Rpad_run_next_node(node);  
//  Rpad_calculate_tree(node);  
}

function Rpad_calculate_tree(node) {
  if (node == null || node.nodeType == 3) return;
  Rpad_calculate_node(node);
  for (var chld = node.firstChild;chld;chld=chld.nextSibling) 
    Rpad_calculate_tree(chld);  
}

function Rpad_run_next_node(node,doit) {
  if (node == null) return;
  if (node.firstChild != null) {
    if (Rpad_calculate_node(node.firstChild,doit)) return;
    Rpad_run_next_node(node.firstChild,doit);
  }
  else { // try siblings and parent's siblings
//alert(node.tagName);
    while (node.nextSibling == null && (node.nodeType != 1 || node.nodeName != 'BODY')) 
      node = node.parentNode;
    if (node.nodeName == 'BODY') // the end of the tree traversal
      _runState = "normal";  
    else { // found a way to keep traversing the tree
      if (Rpad_calculate_node(node.nextSibling,doit)) return;
      Rpad_run_next_node(node.nextSibling,doit);
    }
  }
}
function Rpad_run_next_node(node,doit) { // try at a non-recursive version
  while (node != null) {
    if (node.firstChild != null) {
      if (Rpad_calculate_node(node.firstChild,doit)) return;
      node = node.firstChild;
    }
    else { // try siblings and parent's siblings
      while (node.nextSibling == null && (node.nodeType != 1 || node.nodeName != 'BODY')) 
        node = node.parentNode;
      if (node.nodeName == 'BODY') { // the end of the tree traversal
        _runState = "normal";  
        return;
      }
      else { // found a way to keep traversing the tree
        if (Rpad_calculate_node(node.nextSibling,doit)) return;
        node = node.nextSibling;
      }
    }
  }
}

function inputCommandText(node) {
  var name = node.getAttribute("name");
  if (name == "") return "";
  if (node.type.toLowerCase() == "checkbox")
    if (node.checked)
      return name + " = TRUE";
    else
      return name + " = FALSE";
  else if (node.type.toLowerCase() == "radio" && node.checked) 
    return name + " = '" + node.value + "'";
  else if (node.nodeName.toLowerCase() == "select")
    return name + " = '" + node[node.selectedIndex].text.replace(/'/g,"\\\'") + "'"
//    return "assign('" + name + "', '" + node[node.selectedIndex].text + "')";
  else if (node.value != "")
    if (node.getAttribute("Rpad_type") == "Rstring") {
//        alert(node.value.replace(/'/g,"\\\'"));
        return name + " = '" + node.value.replace(/'/g,"\\\'") + "'"
//      return "assign('" + name + "', '" + node.value + "')";
    }
    else if (node.getAttribute("Rpad_type") == "Rvariable")
      return name + " = " + node.value 
//      return "assign('" + name + "', " + node.value + ")";
  return "";
}

function Rpad_calculate_node(node,doit) { // returns true if a "calculatable" node is found
    var result = "";
    if (node.nodeType != 1) return false;
//alert(node.nodeName+' '+node.nodeType+' '+node.name+' '+node.Rpad_type);
//if (node.nodeName == "INPUT") alert(node.value+'X'+node.nodeType);
      
    var rrun = node.getAttribute("Rpad_run");
    var isReady =  (_runState == "init" && rrun == "init") ||
                   (_runState == "normal"  && (rrun == null || rrun == "normal")) ||
                   (rrun == "all") ||
                   doit;
    if (!isReady) return false;
    if (node.getAttribute("Rpad_type") == "R") {      
      _Rpad_Results = node.getAttribute("Rpad_output");
      R_run_node(getTextContent(node),node);
      return true;
    }
    if (node.getAttribute("Rpad_type") == "shell") {      
      _Rpad_Results = node.getAttribute("Rpad_output");
      shell_run_node(getTextContent(node),node);
      return true;
    }
    if (node.getAttribute("Rpad_type") == "file" && node.getAttribute("name") != null) {
      put_file_node(getTextContent(node),node.getAttribute("name"),node);
      return true;
    }
    if (node.getAttribute("Rpad_type") == "javascript") { 
      // NOTE: javascript only shows the results of the last output
      _Rpad_Results = node.getAttribute("Rpad_output");
      append_results(fixformat(eval(getTextContent(node)).toString(),
                               node.getAttribute("Rpad_output")),
                     node);
    }
    if (node.getAttribute("Rpad_type") == "AsciiMathML") { 
      var n = AMcreateElementXHTML();
      n.appendChild(document.createTextNode("`"+getTextContent(node)+"`"));
      AMinitSymbols();
      AMbody = document.getElementsByTagName("body")[0];
      AMprocessNode(n, false);
      append_results_node(n, node);
    }

    var isRstring = node.getAttribute("Rpad_type") == "Rstring" && 
                    node.getAttribute("name") != null &&
                    node.nodeName != "INPUT" && node.nodeName != "SELECT";
    if (isRstring) {
      _Rpad_Results = "";
      R_run_node(node.getAttribute("name")  + " = '" +
                         getTextContent(node).replace(/\n/g,"\\n") + "'", node);
//      R_run_node("assign('" + node.getAttribute("name")  + "', '" + 
//                        getTextContent(node).replace(/\n/g,"\\n") + "')", node);
      return true;
    }
    var isRvariable = node.getAttribute("Rpad_type") == "Rvariable" && 
                      node.getAttribute("name") != null &&
                      node.nodeName != "INPUT" && node.nodeName != "SELECT";
    if (isRvariable) {
      _Rpad_Results = "none";
      R_run_node(node.getAttribute("name")  + " = " +
                         getTextContent(node).replace(/\n/g,"\\n"), node);
//      R_run_node("assign('" + node.getAttribute("name")  + "', " + 
//                    getTextContent(node).replace(/\n/g,"\\n") + ")", node);
      return true;
    }

    // for forms, process the whole form at once
    var commandtext = "";
    if (node.nodeName == "FORM") 
      for (var elem = 0; elem < node.length; elem++) 
        commandtext = commandtext + inputCommandText(node.elements[elem]) + "\n";

    // process standalone input fields individually
    // assign input fields to language variables
    if ((node.nodeName == "INPUT" || node.nodeName == "SELECT") && 
        node.getAttribute("name") != "" && 
//        node.getAttribute("Rpad_type") != null && 
        hasNoFormParent(node)) 
      commandtext = inputCommandText(node);
    if (commandtext != "") {
//alert(commandtext);
      _Rpad_Results = "none";
      R_run_node(commandtext,node);
      return true;
    }
    var isLoneResult = (node.className == "Rpad_results" && node.previousSibling == null) ||
                       (node.className == "Rpad_results" && 
                        node.previousSibling.className != "Rpad_input" &&
                        node.previousSibling.className != "wrapperForHidden");
//    if (isLoneResult)
//      node.parentNode.removeChild(node);
    return false;
}

function hasNoFormParent(node) {
  while (node.nodeName != "FORM" && node.nodeName != "BODY")
    node = node.parentNode;
  return node.nodeName == "BODY";
}

function Rpad_save() {
  var contents = encodeURIComponent(_Rpad_editor.getHTML().replace(/\r\n/g,"\n"));
//  var contents = encodeURIComponent(_Rpad.innerHTML);
  var filename = location.pathname.split('/'); 
  filename = filename[filename.length-1]; // the last part is the actual file name
  var xmlhttp = myXmlHttp();
  if(Rpad_dir == "") {
    alert("Can't connect to the server");
    return;
  }
  window.status = 'Sending file to the server';
  xmlhttp.open("POST", "server/Rpad_process.pl",true);
  xmlhttp.onreadystatechange=function() {
    if (xmlhttp.readyState==4) {
      if (xmlhttp.responseText.match(/Rpad ERROR:/)) {
        alert(xmlhttp.responseText);
        return;
      }
      var result = fixformat(trim(xmlhttp.responseText));
      window.status = 'Rpad page saved';
    }
  }
  xmlhttp.setRequestHeader("Content-Type", "application/x-www-form-urlencoded"); 
  var s = "ID="+Rpad_dir+"&command=savehtml&filename=" + filename + "&content=" + contents; 
  xmlhttp.send(s); 
}

function Rpad_save_as() {
  var contents = encodeURIComponent(_Rpad_editor.getHTML().replace(/\r\n/g,"\n"));
//alert(contents);
//alert(_Rpad_editor._doc.body.innerHTML);
//  var contents = encodeURIComponent(_Rpad.innerHTML);
  var filename = prompt('SAVE: Enter the file name (without extension)','test') + ".Rpad";
  var xmlhttp = myXmlHttp();
  if(Rpad_dir == "") {
    alert("Can't connect to the server");
    return;
  }
  window.status = 'Sending file to the server';
  xmlhttp.open("POST", "server/Rpad_process.pl",true);
  xmlhttp.onreadystatechange=function() {
    if (xmlhttp.readyState==4) {
      var result = fixformat(trim(xmlhttp.responseText));
      window.status = 'Rpad page saved';
      window.open(filename);
    }
  }
  xmlhttp.setRequestHeader("Content-Type", "application/x-www-form-urlencoded"); 
  var s = "ID="+Rpad_dir+"&command=savehtml&filename=" + filename + "&content=" + contents; 
  xmlhttp.send(s); 
}

function startRpad(){
  Rpad_login();
//  R_login();
  initEditor();
  _runState = "init";  
  setTimeout(function() {
//    document.body.scroll="no";  // not needed: done in htmlarea.css style sheet 
    Rpad_calculate();
  }, 400);  // give it some time
}

function stopRpad(){
  R_logout();
}

