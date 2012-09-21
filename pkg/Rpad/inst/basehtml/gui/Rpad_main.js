/*
    Rpad.js  --  Rpad client-side functionality

    by Tom Short, EPRI, tshort@epri.com

    (c) Copyright 2005 - 2006. by EPRI, Inc.

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


////////////////////////////
// This provides the client-side Rpad functionality. The code is arranged by 
// module function. Basically each module should only depend on functions and 
// objects defined in modules above the given object (and not on modules below).
// The main modules are:
//    dojo.widget.HtmlRpad - the code/results widget that goes on a page
//    rpad - basic code for logging in, sending commands to the server and receiving results
//    rpad.script - registration and processing for various scripts (R, javascript, etc.)
//    rpad calculations - page traversal and calculation control code
//    rpad.utils - these utility functions are defined throughout this 
//                 file, mainly in the order where they're needed.


////////////////////////////
// Load dojo
// 
//
//if (typeof dojo=="undefined") {
//    document.write('<script src="gui/dojo.js" type="text/javascript"></script>');
//}


////////////////////////////
// Rpad configuration options
// 
//

if(dj_undef("rpadConfig")){
    var rpadConfig = {};
}

// default configuration options
var _config = {
    gui: "none",     // options are "none", "standard", "alternate"
    // default options for Rpad widgets:
    rpadType: "R",
    rpadHideSource: false,
    rpadOutput: "normal",   // options are "normal", "none", "javascript"
    rpadRun: "normal",
    onBeforeUnload: false  // ask before leaving a page?
}

if (typeof rpadConfig == "undefined") { 
    rpadConfig = _config; 
} else {
    for (var _option in _config) {
        if (typeof rpadConfig[_option] == "undefined") {
            rpadConfig[_option] = _config[_option];
        }
    }
}
   
////////////////////////////
// The main namespace for Rpad
// 
//

var rpad = {};
   
////////////////////////////
// Utility functions
// Internal utility functions for Rpad
// More are defined throughout this file

rpad.utils = {};
rpad.utils.getTextContent = function(node) {
  if (node.nodeType == 3) { // text node
    if (node.parentNode.nodeName == "PRE") { // really should see if the style has "white-space: pre"
      return node.nodeValue;
    } else {
      return node.nodeValue.replace(/\n/g," ");
    }
  }
  if (node.nodeType == 1) { // element node
    if (node.nodeName == "BR") return "\n";
    if (node.nodeName == "TEXTAREA") return node.value;
    if (node.nodeName == "INPUT") return node.value;
    var text = [];
    for (var chld = node.firstChild;chld;chld=chld.nextSibling) {
        text.push(rpad.utils.getTextContent(chld));
    } return text.join("");
  } return ""; // some other node, won't contain text nodes.
}  

////////////////////////////
// This provides the basic Rpad widget, with input (DOM node rpadInput)
// and output (DOM node rpadResults)
// It's public methods are:
//    initialize(): get it ready
//    toggleVisibility(): hide and unhide source
//    show(): show source
//    hide(): hide source
//    calculate(): calculate it
// The main CSS classes available for modification are:
//    Rpad_input
//    Rpad_output
   
dojo.provide("dojo.widget.Rpad");
//dojo.provide("dojo.widget.HtmlRpad");
   
//dojo.require("dojo.lang.*");
//dojo.require("dojo.io.*");
// dojo.require("dojo.widget.*");
// dojo.require("dojo.event.*");
// dojo.require("dojo.xml.*");
// dojo.require("dojo.xml.Parse");
// dojo.require("dojo.widget.Parse");

dojo.require("dojo.widget.*");
dojo.require("dojo.widget.HtmlWidget");
dojo.require("dojo.event.*");
dojo.require("dojo.html.style");
//dojo.require("dojo.html.selection");
//dojo.require("dojo.debug.console");

rpad.utils.parseNode = function(node) {
// parse the block for other dojo stuff
// bits from David Schontzler 
    var frag = new dojo.xml.Parse().parseElement(node, null, true);
    dojo.widget.getParser().createComponents(frag);
}

rpad.utils.parseNode = function(node) {
// parse the block for other dojo stuff http://dojotoolkit.org/pipermail/dojo-interest/2006-August/014495.html
     dojo.widget.createWidget(node);
}


dojo.widget.defineWidget(
    "dojo.widget.Rpad",
    dojo.widget.HtmlWidget,
{
    templateString:
        '<div dojoAttachPoint="rpadWrapper">                                                        '+
        '    <div dojoAttachPoint="rpadInputWrapper">                                               '+
        '        <span contentEditable = "false">                                                   '+
        '            <textarea dojoAttachPoint="textarea" style="display: none;"></textarea>        '+
        '        </span>                                                                            '+
        '    </div>                                                                                 '+
        '    <div contentEditable = "false">                                                        '+
        '        <div class="Rpad_results" dojoAttachPoint="rpadResults" style="display: none;">    '+
        '            Results                                                                        '+
        '        </div>                                                                             '+
        '    </div>                                                                                 '+
        '</div>                                                                                     ',
                                                                                                   

//    this.templateCssPath = dojo.uri.dojoUri("src/widget/templates/Rpad.css");
   
//    widgetType: "Rpad",
    rpadType: rpadConfig.rpadType,
    rpadOutput: rpadConfig.rpadOutput,
    rpadRun: rpadConfig.rpadRun,
    rpadHideSource: rpadConfig.rpadHideSource,
    visible: true,
   
    tagName: "dojo:rpad",
   
    fillInTemplate: function(args, frag){
         this.rpadInput = frag[this.tagName].nodeRef.cloneNode(true);
         dojo.html.addClass(this.rpadInput, "Rpad_input");
         this.rpadInput.dojotype = null;
         this.rpadInput.removeAttribute("dojotype");
         this.rpadInputWrapper.appendChild(this.rpadInput);
//         rpad.utils.parseNode(this.rpadInput);
   
        if (this.rpadInput.nodeName == "TEXTAREA" || this.rpadInput.nodeName == "INPUT") {
            this.rpadInputWrapper.setAttribute("contentEditable", false);
            this._isNormallyContentEditable = false;
        } else {
            this._isNormallyContentEditable = null;
        }
    },
    initialize: function() {
        this.domNode.rpadType = this.rpadType;
        this.domNode.setAttribute("rpadType", this.rpadType);
        this.domNode.setAttribute("rpadRun", this.rpadRun);
        if (this.rpadHideSource) this.hide();
        this.rpadInput.rpadWidget = this;
        this.domNode.rpadWidget = this;
    },
    show: function() {
        this.rpadInput.style.display = "";
        this.rpadInput.parentNode.contentEditable = this._isNormallyContentEditable; // broken in IE
        this.visible = true;
    },
    hide: function() {
        this.rpadInput.style.display = "none";
        this.rpadInput.parentNode.contentEditable = false;
        this.visible = false;
    },
    toggleVisibility: function() {
        if (this.visible) {
            this.hide();
        } else {
            this.show();
        }
    },
    calculate: function() {
        var rpadResults = this.rpadResults;
        if (this.rpadOutput == "javascript") { // use for JSON or other server-side javascript
            rpadResults = "javascript";
        } else if (this.rpadOutput == "none") {
            rpadResults = null;
        }
        rpad.send(this.rpadType, rpad.utils.getTextContent(this.rpadInput), 
                  rpadResults, this.rpadInput);
    }, 
    onReceive: function() { // stub available for attaching
    },
    getHtml: function() { // get a HTML representation of the widget
        var result = "<" + this.rpadInput.tagName.toLowerCase() + ' dojoType="Rpad"';
        if (this.rpadType != rpadConfig.rpadType)
            result += ' rpadType="' + this.rpadType + '"';
        if (this.rpadOutput != rpadConfig.rpadOutput)
            result += ' rpadOutput="' + this.rpadOutput + '"';
        if (this.rpadRun != rpadConfig.rpadRun)
            result += ' rpadRun="' + this.rpadRun + '"';
        if (this.rpadHideSource != rpadConfig.rpadHideSource)
            result += ' rpadHideSource="' + this.rpadHideSource + '"';
        result += ">" + rpad.utils.getTextContent(this.rpadInput) + 
                  "</" + this.rpadInput.tagName.toLowerCase() + ">";
        return result;
    }
})  
   
// dojo.inherits(dojo.widget.HtmlRpad, dojo.widget.HtmlWidget);
// dojo.widget.tags.addParseTreeHandler("dojo:rpad");

////////////////////////////
// NOT USED!!!
// This provides a basic Rpad form widget. It's a small wrapper around
// input's, textarea's, and selects.
// It's public methods are:
//    calculate(): calculate it
//    onDoneCalc(): calculations done
   

dojo.widget.RpadForm = function(){
    dojo.widget.Widget.call(this);
   
    this.widgetType = "RpadForm";
    this.tagName = "dojo:rpadform";
   
    this.calculate = function() {
        rpad.send(this.rpadType, rpad.utils.getTextContent(this.domNode), 
                  this.rpadResults, this.rpadInput);
    } 
    this.onDoneCalc = function() { // stub available for attaching
    }
}  
   
// dojo.inherits(dojo.widget.RpadForm, dojo.widget.HtmlWidget);
// dojo.widget.tags.addParseTreeHandler("dojo:rpadform");
   
////////////////////////////
// This provides the base Rpad class (it's not R specific)
// It's public methods are:
//    send(rpadType, commands, rpadResults, node): send for processing
//    receive(rpadResults, data, startNode): receive processed results
//    updateResults(rpadResults, data, rpadOutputStyle): update results received
// Note: More methods are defined below.
   
rpad.dir = "";
rpad._runState = "init";  
rpad._activeWindow = window; 
   
rpad.send = function(rpadType, commands, rpadResults, rpadInput, jsOnDone){
// Send "commands" for processing. 
// "commands" are script text of type "rpadType": "R", "file", ...
// Put results received in the DOM node "rpadResults".
// Send along the originating DOM node "rpadInput" in case it's needed.
//dojo.debug('SEND:'+commands);
    rpad.script.run(rpadType, commands, rpadResults, rpadInput, jsOnDone)
//    dojo.debug('SEND: '+commands);
}  
   
rpad.sendThenRun = function(rpadType, commands, jsOnDone){
// NOT WORKING YET
// Send "commands" for processing. 
// "commands" are script text of type "rpadType": "R", "file", ...
// Ignore results.
// Eval the "jsOnDone" when results received.
//dojo.debug('SEND:'+commands+rpadInput.nodeName);
  rpad.script.run(rpadType, commands, rpadResults, rpadInput)
}  
   
rpad.receive = function(rpadResults, data, rpadInput, jsOnDone) {
// Receive the processed commands "data" (text/html).
// Put results received in the DOM node "rpadResults".
// Send along the originating DOM node "rpadInput" in case it's needed.
// if rpadResults == "javascript", exec the javascript instead of inserting results in the DOM
//dj_debug("REC:"+rpadInput.nodeName+data);
//dj_debug("REC:"+data);
//    dojo.log.debug('REC:\n'+data);
    if (rpadInput.onReceive) rpadInput.onReceive(); // fire an event inputs can attach to
    if (rpadInput.rpadWidget) {
        rpadInput.rpadWidget.onReceive(); // fire an event Rpad widgets can attach to
        var rpadOutputStyle = rpadInput.rpadWidget.rpadOutput;
    }
    if (rpadResults == "javascript") {
        dj_eval(data);
    } else if (rpadResults != null) {
        rpadResults.style.display = "";
        rpad.updateResults(rpadResults, data, rpadOutputStyle);       
    }
    if (jsOnDone) 
        setTimeout(function() {dj_eval(jsOnDone)}, 50); // IE7 can die without the timeout
    if (rpad._doKeepGoing) { // this would be better separated, but connect was acting funny
        rpad.calculateNextNode(rpadInput);
    }
}  
   
rpad.updateResults = function(rpadResults, data, rpadOutputStyle) {
// Insert the processed commands "data" (text/html)
// in the DOM node "rpadResults".
// Decode the text/html and handle the switching between
// html mode and text mode (done with <htmlon/> and <htmloff/>.
    function fixhtmlencodings(str) {  
      str = str.replace(/&/g,"&amp;");  
      str = str.replace(/</g,"&lt;");  
      str = str.replace(/>/g,"&gt;");  
      str = str.replace(/\n/g,"\n<BR/>");
      str = str.replace(/ /g,"&#160;"); // replace spaces with nonbreaking spaces
      return(str);
    }
    
    function fixformat(str, type) {
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
    if (dojo.render.html.ie) {
        rpadResults.innerHTML = fixformat(data, rpadOutputStyle);  
    } else { // move rpadResults out of the DOM and back in to get FF to execute embedded scripts
        var marker = document.createElement("DIV");
        dojo.dom.insertAfter(marker, rpadResults);
        var node = document.createElement("DIV");
        node.appendChild(rpadResults);
        rpadResults.innerHTML = fixformat(data, rpadOutputStyle);  
        marker.parentNode.replaceChild(rpadResults, marker);
    }    
    rpad.utils.parseNode(rpadResults);
}  
   
////////////////////////////
// Rpad login: send a login command to the server and update the 
// appropriate working directory.
// Hook into this to add other login handlers.
   
rpad.login = function() {
    dojo.io.bind({
        url: "server/Rpad_process.pl?command=login",
        load: rpad.afterLogin,
        error: function(type, error){ dojo.debug(error); },
        mimetype: "text/html"
    });
}  
rpad.afterLogin = function(type, data, evt){
    rpad.dir = data;
}

   
////////////////////////////
// This provides script registration/processing
// It's public methods are:
//    register(rpadType, run): register the given js function "run" of type "rpadType"
//    unregister(rpadType): unregister the script "rpadType"
//    run(rpadType, commands, rpadResults, node): run the given script
//    serverScript(url, getContent): return a js function that will send data to a server with the appropriate url
//                          "getContent" is a function with arguments (commands, node) 
   
rpad.script = {};
rpad.script._scripts = [];
rpad.script.registered = function() {
    return rpad.script._scripts;
}  
rpad.script.register = function(rpadType, run) {
    rpad.script._scripts[rpadType] = run;
}  
rpad.script.unregister = function(rpadType) {
    rpad.script._scripts[rpadType] = null;
}  
rpad.script.run = function(rpadType, commands, rpadResults, rpadInput, jsOnDone) {
    rpad.script._scripts[rpadType](commands, rpadResults, rpadInput, jsOnDone); 
}  
rpad.script.serverScript = function(url, getContent) {
// Returns a js function that will send data to a server 
// with the appropriate url.
// "getContent" is a function with arguments (commands, rpadResults, node) that returns
// the content to be passed to the server.
    return function(commands, rpadResults, rpadInput, jsOnDone) {    
        dojo.io.bind({
            url: url,
            content: getContent(commands, rpadInput),
            method: "POST",
            load: function(type, data, evt){
                 rpad.receive(rpadResults, data, rpadInput, jsOnDone)
            },
//            error: function(type, error){ dojo.debug(error); },
            error: function(type, error){ props(error); dojo.debug(error.number); dojo.debug(error.message);},
            mimetype: "text/html"
        });
    }
}  

// This provides a stub for form processing:
   
rpad.processForm = function(nodes) {};
  
   
////////////////////////////
// Rpad calculation routines, providing the following methods that the user
// could call:
//   - calculatePage(): Calculates the whole page.
//   - calculateNext(node): Calculates the first Rpad widget after or under the given DOM node 
// The following are used internally but may be "connected to" to change functionality:
//   - calculateNextNode(node, doit): Calculates the first Rpad widget after the given DOM node 
//   - calculateNode(node, doit): Calculates the given DOM node (Rpad widget or form element)
// The page also keeps track of the internal state of the page with the following:
//   - _doKeepGoing: If true, keeps traversing the DOM tree and calculating appropriate nodes
//   - _runState: "init" initially and "normal" after the first pass through the page
   
   
rpad.calculatePage = function(){
//dojo.debug('calcPage');
  rpad.calculateTree(rpad.base);  
}  
rpad.calculateTree = function(node) {
//dojo.debug("CALCTREE: "+node.tagName+node.parentNode.tagName+node.parentNode.parentNode.tagName);
  rpad._doKeepGoing = true;  // turn on page traversal
  rpad._serverTries = 0;
  rpad._startingNode = node;  // remember the starting point
  rpad.calculateNextNode(node);  
}  
rpad.afterPageCalculation = function() {
//dojo.debug("afterPageCalculation");
    rpad._doKeepGoing = false;  // turn off page traversal
    rpad._runState = "normal";  
}; 
rpad.pageError = function() {dojo.debug('page error')};

rpad.calculateNextNode = function(node, doit) { // non-recursive
// Calculate the next Rpad element after or below the DOM node "node".
// Traverses the DOM tree starting at the DOM node "node".
// Keeps going until it finds and calculates a "calculatable" DOM node.
// "doit" means calculate it regardless of rpad._runState and the status of node.rpadRun
//  try {
    while (typeof(node) != "undefined" && node != null) {
      if (node.firstChild != null) { // try children
        if (rpad.calculateNode(node.firstChild, doit)) {
          return;
        }
        node = node.firstChild;
      }
      else { // try siblings and parent's siblings
      while (node.nextSibling == null && (node.nodeType != 1 || node.nodeName != 'BODY') && node != rpad._startingNode) { 
          node = node.parentNode;
        }   
        if (node.nodeName == 'BODY' || node == rpad._startingNode) { // the end of the tree traversal
          rpad.afterPageCalculation();
          return;
        } else { // found a way to keep traversing the tree
          if (rpad.calculateNode(node.nextSibling,doit)) {
            return;
          }
          node = node.nextSibling;
        }
      }
    }
//  } catch(err) {
//    rpad.pageError();
//  }
  rpad.afterPageCalculation();
}  
   
rpad.utils.getWidget = function(node) {
// returns the dojo widget associated with the given DOM "node"
// Dylan Schiemann
// http://dojotoolkit.org/pipermail/dojo-interest/2005-September/000913.html
    var unaryFunc = function(x){
        if(x.domNode == node){
            return true;
        }
    }
    return dojo.widget.getWidgetsByFilter(unaryFunc)[0];
}  
rpad.utils.getRpadWidget = function(node) {
// returns the dojo widget associated with the given DOM "node"
// Dylan Schiemann
// http://dojotoolkit.org/pipermail/dojo-interest/2005-September/000913.html
    var unaryFunc = function(x){
        if(x.domNode == node && x.widgetType == "Rpad"){
            return true;
        }
    }
    return dojo.widget.getWidgetsByFilter(unaryFunc)[0];
}  
rpad._hasNoFormParent = function(node) {
// returns true if the DOM "node" does not have a FORM-element parent
    while (node.nodeName != "FORM" && node.nodeName != "BODY")
        node = node.parentNode;
    return node.nodeName == "BODY";
}  
   
rpad.calculateNode = function(node, doit) { // returns true if a "calculatable" node is found
// Calculates the DOM node "node".
// "doit" means calculate it regardless of rpad._runState and the status of node.rpadRun
// Returns true if it found a node to calculate, false if not.
    if (node.nodeType != 1) return false;
      
    var rrun = node.getAttribute("rpadRun");
    var isReady =  (rpad._runState == "init" && rrun == "init") ||
                   (rpad._runState == "normal"  && (rrun == null || rrun == "" || rrun == "normal")) ||
                   (rrun == "all") ||
                   doit;
//    if (!isReady) return false;
    if (typeof( node.rpadType ) != "undefined" && isReady) {      
//dojo.debug("CALCNODE: "+node.tagName+node.parentNode.tagName+node.parentNode.parentNode.tagName+node.rpadType);
        var widget = rpad.utils.getWidget(node);
        if (typeof(widget) != "undefined" && widget.calculate != null) {
            widget.calculate();
            return true;
        }
    }
    // for forms, process the whole form at once
    if (node.nodeName == "FORM") {
        rpad.processForm(node.elements);
        return true;
    }
    // process standalone input fields individually
    if ((node.nodeName == "INPUT" || node.nodeName == "SELECT") && 
        node.getAttribute("name") != "" && 
        rpad._hasNoFormParent(node)) {
            rpad.processForm(node);
            return true;
    }
    return false;
}  
   
// Mapping for backwards compatibility
var Rpad_calculate = rpad.calculatePage;
   
// make Rpad keep calculating nodes after each node is processed
// had to disable it and add it directly to rpad.receive because it was acting funny on Mozilla
//dojo.event.connect(rpad, "receive", 
//    function (rpadResults, data, rpadInput) {
//        if (rpad._doKeepGoing)
//            rpad.calculateNextNode(rpadInput);
//    })
   
rpad.calculateNext = function(node) { 
// Calculate the next Rpad element after or below the DOM node "node".
// Make sure it doesn't keep going, and 
// calculate even if it's an init block.
   rpad._doKeepGoing = false;  // turn off page traversal
   rpad.calculateNextNode(node, true);
}  

////////////////////////////
// Key handling
//   just the basics here - F9 for calculatePage
   
rpad.onKey = function(e) {
    if (e.key == e.KEY_F9) {
        rpad.calculatePage();
        dojo.event.browser.stopEvent(e);
    }
    if (e.key == e.KEY_F11 && e.ctrlKey) { 
        // toggle debug:
        djConfig.isDebug = !djConfig.isDebug;
    }
    if (e.key == e.KEY_F12 && e.ctrlKey) {
        if (rpad.currentGui == "none") {
            rpad.gui2.init();
        } else {
            rpad.gui2.disable();
        }
        dojo.event.browser.stopEvent(e);
    }
}

////////////////////////////
// Loading... display
// 
rpad.addLoadingDotDot = function() {
    rpad.loadingMessageNode = document.createElement('div');
    rpad.loadingMessageNode.innerHTML = "<img src='gui/wait.gif'> Page loading...";
    rpad.loadingMessageNode.id = "rpadPageLoadingMessage"
    document.body.appendChild(this.loadingMessageNode);
}   
rpad.removeLoadingDotDot = function() {
    dojo.dom.removeNode(rpad.loadingMessageNode);
}   


////////////////////////////
// R-specific initialization
//   - install R scripts
//   - install an R form handler
//   - install an R login handler (that starts R)
   
rpad.script.register("R", 
    rpad.script.serverScript("server/R_process.pl",
        function(commands, rpadInput) { 
            return {ID: rpad.dir, 
                    command: "R_commands", 
                    R_commands: commands}
        }));
   
//rpad.script.register("Rstring", 
//    rpad.script.serverScript(
//        function(commands, rpadInput) { 
//            return "server/R_process.pl?" + 
//                   "ID="+rpad.dir+
//                   "&command=R_commands&R_commands="+
//                   rpadInput.name+"='"
//                   encodeURIComponent(commands) + "'";}))
//rpad.script.register("Rvariable", 
//    rpad.script.serverScript(
//        function(commands, rpadInput) { 
//            return "server/R_process.pl?" + 
//                   "ID="+rpad.dir+
//                   "&command=R_commands&R_commands="+
//                   rpadInput.name+"="
//                   encodeURIComponent(commands);}))
   
rpad.processRForm = function(nodes) {
// Sweep through the elements of the DOM form "node" and
// generate a list of assignment commands for R.
// Then send the commands to R.
    var commands = "";
    if (!nodes.length || (nodes.nodeName && nodes.nodeName.toLowerCase() == "select")) 
        nodes = [nodes]; // make individual nodes into an array
    for (var i=0; i < nodes.length; i++) {
        var node = nodes[i];
        var name = node.getAttribute("name");
//dojo.debug(name+"X"+node.value);
        if (name == "") continue;
        var command = "";
        if (node.type == "checkbox") {
          if (node.checked)
            command = name + " = TRUE";
          else
            command = name + " = FALSE";
        } else if (node.type == "radio") {
          if (node.checked) 
            command = name + " = '" + node.value + "'";
        } else if (node.nodeName.toLowerCase() == "select" && node.selectedIndex >= 0)
          command = name + " = '" + node[node.selectedIndex].text.replace(/'/g,"\\\'") + "'"
        else if (node.type == "text" || node.type == "hidden")
          if (node.getAttribute("rpadType") == "Rvariable" && node.value != "") {
            command = name + " = " + node.value; 
          }
          else
            command = name + " = '" + node.value.replace(/'/g,"\\\'") + "'";
        commands = commands + command + "\n";
    }
    if (commands != "") {
        rpad.send("R", commands, null, nodes[0]);
    }
}  
// Connect the R form handler
dojo.event.connect(rpad, "processForm", rpad, "processRForm");
  
// Connect the R login handler
rpad.loginR = function() {
    dojo.io.bind({
        url: "server/R_process.pl?ID="+rpad.dir+"&command=login",
        load: rpad.afterLoginR,
        error: function(type, error){ dojo.debug(error); },
        mimetype: "text/html",
        timeoutSeconds: 3,         // added to make R "load-and-go" to make it work better with CGI
        timeout: rpad.afterLoginR  // just assume it worked okay
    });
}  
rpad.afterLoginR = function() {
    rpad._isRReady = true;
    rpad.calcIfReadyR();
}
rpad.calcIfReadyR = function() {
    if (rpad.pageIsLoaded && rpad.dir != "" && rpad._runState == "init" && rpad._isRReady)
        rpad.calculatePage();
}

rpad.logoutR = function() {
    dojo.io.bind({
        url: "server/R_process.pl?ID="+rpad.dir+"&command=logout"
    });
}  
dojo.event.connect("before", window, "onunload", rpad, "logoutR");

//dojo.event.connect(rpad, "afterLogin", rpad, "afterLoginR"); // attempt to fix a cludgy bug
// we need both of these to happen, but we don't know which will happen first:
dojo.event.connect(rpad, "afterLogin", rpad, "loginR");

//dojo.event.connect(dojo, "loaded", rpad, "addLoadingDotDot");
//dojo.event.connect(rpad, "calculatePage", rpad, "removeLoadingDotDot");
//dojo.event.connect(rpad, "afterPageCalculation", rpad, "removeLoadingDotDot");

////////////////////////////
// Initialization
//   * log in right away   
//   * set up keyboard and login events
//   * The order logic is: 
//          rpad.login - rpad.afterLogin - rpad.loginR - rpad.afterLoginR - calculatePage
//                                                                               | 
//          dojo.loaded - rpad.init - rpad.afterLoginR ---------------------------
//
//          
   
rpad.pageIsLoaded = false;  // used to determine if the page is ready for a calculation
   
rpad.login(); // login right away -- don't wait for the page to finish loading

rpad.init = function(e) {
    rpad.currentGui = "none";
    rpad.pageIsLoaded = true;
//    rpad.calcIfReadyR();
    rpad.base = document.body; // reference point for all calculations
    //rpad.login(e);
    dojo.event.connect(document, "onkey", this, "onKey");
//     dojo.event.browser.addListener(document, "onKey", dojo.lang.hitch(this, this.onKey));
//     var events = ["onkeydown", "onkeypress"];
//     for (var i=0; i < events.length; i++) {
//         dojo.event.connect(document, events[i], rpad, "handleEvents");
//     }
    if (rpadConfig.onBeforeUnload) window.onbeforeunload = function() {return "";};
}  
//dojo.event.connect(dojo, "loaded", function() {dj_debug("loaded")})
dojo.event.connect(dojo, "loaded", rpad, "init")
dojo.event.connect(dojo, "loaded", rpad, "calcIfReadyR"); // this is out of place, but it won't work if above, either IE or Moz won't autoload.

   
////////////////////////////
// Other script registrations
//   - shell, file, javascript
   
rpad.script.register("shell", 
    rpad.script.serverScript("server/shell_process.pl",
        function(commands, rpadInput) { 
            return {ID: rpad.dir, 
                    shell_commands: commands}
        }));
rpad.script.register("file", 
    rpad.script.serverScript("server/Rpad_process.pl",
        function(commands, rpadInput) { 
            return {ID: rpad.dir, 
                    command: "savefile",
                    filename: rpadInput.getAttribute("filename"),
                    content: commands}
        }));
rpad.script.register("javascript", 
    function(commands, rpadResults, rpadInput) { 
        rpad.receive(rpadResults,
                          eval(commands).toString(), rpadInput);
    });
      

////////////////////////////
// Throbber to show a page (or tree) calculation in progress
// 
   
rpad.addPageThrobber = function() {
    rpad.throbber = document.createElement('div');
    rpad.throbber.innerHTML = "<img src='gui/wait.gif'>";
    rpad.throbber.id = "rpadPageThrobber"
    document.body.appendChild(rpad.throbber);
}   
rpad.hidePageThrobber = function() {
    if (rpad.throbber)
        rpad.throbber.style.display = "none";
}   
rpad.showPageThrobber = function() {
    if (rpad.throbber)
        rpad.throbber.style.display = "";
}

dojo.event.connect(dojo, "loaded", rpad, "addPageThrobber");
dojo.event.connect(rpad, "calculatePage", rpad, "showPageThrobber");
dojo.event.connect(rpad, "calculateTree", rpad, "showPageThrobber");
dojo.event.connect(rpad, "afterPageCalculation", rpad, "hidePageThrobber");

////////////////////////////
// A drag and drop widget - each drop zone creates a list in R containing the 
// contents of the dropped objects.

dojo.require("dojo.dnd.*");
   
dojo.widget.defineWidget(
    "dojo.widget.RpadDropZone",
    dojo.widget.HtmlWidget,
{
    isContainer: true,
    listName: "browserList",
    dropsAllowed: "*",
    tagName: "dojo:rpaddropzone",
    onDrop: function(e) {},
    afterDrop: function(e) {},
    rpadType: "R",
    rpadRun: rpadConfig.rpadRun,
    
    initialize: function() {
        this.dropTarget = new dojo.dnd.HtmlDropTarget(this.domNode, this.dropsAllowed.split(','));
        dojo.event.connect(this.dropTarget, "onDrop", this, "onDrop");
        dojo.event.connect(this, "onDrop", this, "afterDrop");
        this.domNode.rpadType = this.rpadType;
        this.domNode.setAttribute("rpadType", this.rpadType);
        this.domNode.setAttribute("rpadRun", this.rpadRun);
    },
    calculate: function() {
        var command = this.listName + "= list("; 
        var isListStarted = false;
        for (var i = 0; i < this.domNode.childNodes.length; i++) {
          var chld = this.domNode.childNodes[i];
          if (chld.nodeType != 1) continue;
          var widget = rpad.utils.getWidget(chld);
          if (widget.varValue != null) {
            if (isListStarted) 
              command += ",";
            if (widget.varName != null)
              command += widget.varName + "=";
            command += "'" + widget.varValue.replace(/'/g,"\\\'") + "'"; //"
            isListStarted = true;
          }
        }
        command += ")";
        rpad.send(this.rpadType, command, null, this.domNode);
    }
})  

dojo.widget.defineWidget(
    "dojo.widget.RpadDropZoneSingle",
    dojo.widget.RpadDropZone,
{
    baseZone: "",
    onDrop: function(e) {
        // move all but the dragSource to the baseZone
        var srcNode = e.dropTarget.parentNode;
        var destNode = dojo.byId(this.baseZone);
        while (srcNode.hasChildNodes() && srcNode.firstChild != e.dragSource.domNode) {
            destNode.appendChild(srcNode.firstChild);
        }
        while (e.dragSource.domNode.nextSibling) {
            destNode.appendChild(e.dragSource.domNode.nextSibling);
        }
        rpad.calculatePage();
    }
});

dojo.widget.defineWidget(
    "dojo.widget.RpadDragVariable",
    dojo.widget.HtmlWidget,
{
    varName: null,
    varValue: "varvalue",
    dragName: "drag1",
    tagName: "dojo:rpaddragvariable",
    dropsAllowed: "",
    dropTarget: null,
    initialize: function() {
        this.dragSource = new dojo.dnd.HtmlDragSource(this.domNode, this.dragName);
        this.domNode.innerHTML = this.varValue + " ";
        if (this.dropsAllowed != "") {
          this.dropTarget = new dojo.dnd.HtmlDropTarget(this.domNode, this.dropsAllowed.split(','));
          if (this.onDrop != null)
            this.dropTarget.onDrop = this.onDrop;
        }
    }
}); 

R = function(R_commands){
// Send "R_commands" for processing. Return results as a string
    djConfig.isDebug = true;
    rpad.script.run("R", R_commands, null, document, "djConfig.isDebug = false");
}
   
