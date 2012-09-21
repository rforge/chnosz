guiCommand <- function(funname=NULL, base=NULL, guiHeader=funname, useFunctionDefaults=FALSE,
                       useFunctionParameters=FALSE, ...) {
# Returns a list of class "guiCommand"
# that holds the information needed by HTML.guiCommand to generate a GUI.
# It inherits the definition from "base" if specified.
# It uses defaults from the formal function definition if 'useFunctionDefaults', and
# it shows all function parameters if 'useFunctionParameters'.
  if (is.null(base)) {
    base <- list()
    class(base) <- "guiCommand"
    base$params <- list()
  }
  if (!is.null(funname)) base$funname <- funname
  if (!is.null(guiHeader)) base$guiHeader <- guiHeader
  formls = formals(funname)
  formls = formls[names(formls) != "..."]
  if (useFunctionParameters) {
    for (i in seq(along=formls))
      base$params[[names(formls[i])]] <-
        guiInput(names(formls[i]),
                 default = if (useFunctionDefaults) formls[[i]],
                 type = ifelse(is.character(formls[[i]]),"Rstring","Rvariable"),
                 guiType = if(!is.na(formls[[i]]) && is.logical(formls[[i]])) "guiLogical" else "guiInput"
                 )
  }
      
  params <- list(...)
  for (i in seq(along=params)) {
    base$params[[names(params[i])]] <- params[[i]]
    if (useFunctionDefaults) {
      if (is.null(base$params[[names(params[i])]]$default) &&
          names(params[i]) %in% formls)
        base$params[[names(params[i])]]$default <- formls[[names(params[i])]]
    }
  }
  base
}

guiCommandFun <- function(funname, base=NULL, guiHeader=funname, useFunctionDefaults=FALSE,
                          useFunctionParameters=FALSE, ...) {
# Returns a function that when evaluated returns a list of class "guiCommand"
# that holds the information needed by HTML.guiCommand to generate a GUI.
# The ... parameters to the resultant function are taken as defaults.
  gC = guiCommand(funname=funname, base=base, guiHeader=guiHeader,
    useFunctionDefaults=useFunctionDefaults, useFunctionParameters=useFunctionParameters, ...)
  function(...) {
    params <- list(...)
    for (i in seq(along=params))
      if (names(params[i]) %in% names(gC$params))
        gC$params[[names(params[i])]]$default <- params[[i]]
    gC
  }
}

guiInput <- function(name=NULL, default=NULL, choices=NULL, mapping=NULL, special=NULL, type = "Rvariable", guiType = "guiInput", ...) {
  result <- list(name=name,default=default,choices=choices,mapping=mapping,special=special,type=type) # probably a shortcut for this
  class(result) <- guiType
  result
}
guiInput <- function(name = "", type = "Rvariable", guiType = "guiInput", ...) {
  result <- list(name = name, ..., type = type)
  class(result) <- guiType
  result
}
HTML.guiInput <- function(p, varname) {
  cat("<td>",p$name,"</td><td>")
  print(HTMLinput(id=varname, name=varname, value= p$default, rpad_type=p$type,
                  style='width:160px', contenteditablewrapper=FALSE))
  cat("</td>\n")
}
HTML.guiCombobox <- function(p, varname) {
  cat("<td>",p$name,"</td><td>")
  if (is.null(p$mapping)) choices <- p$choices
  else choices <- doMap(p$choices,p$mapping)
  cat("<span style='position:relative'>")
  print(HTMLinput(id=varname, name=varname, value= p$default, rpad_type=p$type, contenteditablewrapper=FALSE,
                  style="0px;left:0px;z-index:1000;width:160px;"))
  print(HTMLselect(name="dummy", text=choices, contenteditablewrapper=FALSE, optionvalue=p$optionvalue,
                   style="position:absolute;top:0px;left:0px;z-index:0;width:175px;clip:rect(0px,175px, 22px, 160px);",
                   onChange="javascript:var el=this; while(el.nodeName!=\"INPUT\") el = el.previousSibling;el.value = this.options[this.selectedIndex].value;"))
  cat("</span>")
  cat("</td>\n")
}
HTML.guiLogical <- function(p, varname) {
  Html(guiInput(p$name, choices = c(TRUE, FALSE), guiType="guiCombobox"), varname)
}
HTML.guiPickPch <- function(p, varname) {
  cat("<td>",p$name,"</td><td>")
#  cat('<script type="text/javascript">function updateInput(node, val) {')
#  cat(' node.previousSibling.previousSibling.value = val; node.style.display="none";')
#  cat('} </script>\n')
  print(HTMLinput(id=varname, name=varname, value= p$default, rpad_type=p$type,
            style='width:160px', contenteditablewrapper=FALSE))
  cat('<input id="button1" type="button" onclick="javascript:this.nextSibling.style.display=\'\'\n" value="..." />')
  cat('<span  style="display:none; position:relative; background-color:white;"><div style="background-color:white; width:190px; height: 190px; position:absolute;">')
  cat(paste('<img src="js/pch',1:25,'.png" onclick="javascript:updateInput(this.parentNode.parentNode,',1:25,');">\n',sep="",collapse="\n"))
  cat('</div></span>')
  cat("</td>\n")
}

HTML.guiCommand <- function(cmd,varname=make.names(deparse(substitute(cmd))),
                            guiHeader = cmd$guiHeader,
                            allowDuplication = FALSE, expanded = TRUE) {
# Inputs: 
# gC$params$x
# gC$params$y
# gC$params$type
#  
# Creates the following in HTML input boxes: 
# cmd$x
# cmd$y
# cmd$type
#
  gC <- cmd
  params <- gC$params
  cat("<div contenteditable='false'>")
  cat("<a href='\#' onclick='javascript:var el=this; while(el.nodeName!=\"DIV\") el = el.nextSibling; if(el.style.display==\"none\") el.style.display=\"\"; else el.style.display=\"none\";return false;'>",
      guiHeader, "</a><div ", ifelse(expanded,"","style='display:none'"), ">")
  cat("<form>")
# write the hidden values
  # varname <- list()
  print(HTMLinput(name=varname, value="list()", type="hidden", rpad_type="Rvariable"))
  # attr(varname, 'Rcmd') <- capture.output(dput(cmd))
  print(HTMLinput(name=paste("attr(", varname,", \"Rcmd\")"),
                  value=paste(capture.output(dput(cmd)), collapse=" "),
                  type="hidden", rpad_type="Rvariable"))
  idx <- paste(varname,"idx",sep="")
  # idx <- 0
  print(HTMLinput(name=idx, value=0, type="hidden", rpad_type="Rvariable"))
  cat("<span name='cmd_wrapper'><table>")
  # idx <- idx + 1
  print(HTMLinput(name=idx, value=paste(idx,"+ 1"), type="hidden", rpad_type="Rvariable"))
  # varname[[idx]] = list()
  print(HTMLinput(name=paste(varname,'[[',idx,']]',sep=""),
                  value="list()", type="hidden", rpad_type="Rvariable"))
  for (i in seq(along=gC$params)) {
    p <- gC$params[[i]]
    name <- p
    fullname <- paste(varname,'[[',idx,']]$',names(gC$params[i]),sep="")
    cat("<tr>\n")
    Html(p, fullname)
    cat("</tr>\n")
  }
  cat("<td></td>")
  cat("</table>")
  if (allowDuplication) cat("&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp<a href='\#' onClick='javascript:duplicateGui(this);return false;'>duplicate", gC$guiHeader, "</a>")
  cat("</span>")
  cat("</form></div></div>")
}  
doMap <- function(vals,mapping) {
  ifelse(is.na(mapping[vals]), vals, mapping[vals])
}
doInverseMap <- function(val,mapping) {
  ifelse(val %in% mapping,names(mapping[mapping==val]),val)
}

runGuiCommand <- function(cmd,
                          fun = if(!is.null(attr(cmd,'Rcmd')$funname)) get(attr(cmd,'Rcmd')$funname),
                          ...) {
  gC <- attr(cmd,'Rcmd')
  params <- cmd
  for (i in seq(along=gC$params)) {  # invert any mappings
    paramName <- names(gC$params[i])
    mapping <- gC$params[[i]]$mapping
    if (!is.null(mapping))
      for (j in seq(along=params))
        if (!is.null(params[[j]][[paramName]])) 
          params[[j]][[paramName]] <- doInverseMap(params[[j]][[paramName]],mapping)
  }
  for (j in seq(along=params)) {
    do.call("fun", c(params[[j]], ...))
  }
  #  eval(expression(do.call(gC$funname,params)),-2) # evaluate in the routine that calls runGuiCommand
}  



#runGuiCommand(p)
#HTMLon()
#showgraph()


makePchImageMap <- function() {
  require(grid)
  w <- h <- 200
  plot.new()
  png('js/pch.png', width=w*4, height=h*4, pointsize = 12*4) # scale up by a factor of 4
#  bitmap('pch.png', width=w/120, height=h/120, res=120, type="pngalpha", pointsize=12, family="Helvetica") # this didn't look good

  Pex <- 1.5 ## good for both .Device=="postscript" and "x11"
  # stuff taken from ?points
  ipch <- 1:(np <- 25)
  k <- floor(sqrt(np))
  iy <- rev((ipch-1) %/% k + .5) / k
  ix <- rev((k-1)-(ipch-1) %% k + .5) / k
  pch <- as.list(ipch)
  print(ix)
  print(iy)
  for(i in 1:np) {
       pc <- pch[[i]]
       grid.points(unit(ix[i] + .04,"npc"), unit(iy[i],"npc"), pch = pc, size = unit(Pex/30,"npc"), gp = gpar(col = "red", fill = "yellow", lwd=4), default.units="npc")
       ## red symbols with a yellow interior (where available)
       grid.text(pc, ix[i] - .04, iy[i])
     }
  dev.off()
  y1 <- array(rep(seq(1,by=floor(w/k), to=w),k),c(k,k))
  y2 <- y1 + floor(w/k) - 1
  x1 <- t(y1)
  x2 <- t(y2)
  val <- 1:25
  cat('\n<map name="pchmap">\n')
  cat(paste('<area shape="rect" coords="',x1,",",y1,",",x2,",",y2,'" href="javascript:updateInput(\'', val, '\')">\n',collapse="",sep=""))
  cat('</map>\n')
  
}
makePchImages <- function() {
  require(grid)
#  bitmap('pch.png', width=w/120, height=h/120, res=120, type="pngalpha", pointsize=12, family="Helvetica") # this didn't look good

  Pex <- 1.5 ## good for both .Device=="postscript" and "x11"
  # stuff taken from ?points
  w <- h <- 30
  for(i in 1:25) {
    png(paste('js/pch',i,'.png',sep=""), width=w*4, height=h*4, pointsize = 12*4) # scale up by a factor of 4
    par(mar=c(0,0,0,0))
    plot(.75, 0, pch=i, axes=F, xlab="", ylab="", xlim=c(0,1), cex=2, lwd=3, bg="yellow")
    text(.5, 0, labels=i, adj=1)
    dev.off()
  }
}


