\name{json}
\alias{json}
                   
\title{ JSON (JavaScript Object Notation) generation}

\description{
  Write out a JSON structure of an R object.
}

\usage{
  json(x,...)
}

\arguments{
  \item{x}{ an object to convert to JSON.}
  \item{...}{ unused for now.}
}
\details{
JSON is Javascript Object Notation, a data-exchange format especially
suited for passing to Javascript. Mapping from S to Javascript is as
follows:

list (including data frames and other list objects) -> object
character -> string
array  -> array
array with names -> object
matrix -> 2-D array
numeric -> number
NULL -> null
TRUE -> true
FALSE -> false
NA -> NaN

Other types are converted to character using \code{as.character} and
output as character arrays or objects. Handling of NA's is problematic
as Javascript doesn't have NA's. Strings are escaped. Currently,
dimensioning in higher dimension arrays is lost. Object attributes are
also ignored.
}

\value{
  A character string of class json. With automatic printing, the string is sent to the
  output (with \code{cat}).
}

\author{ Tom Short, EPRI, (\email{tshort@epri.com}) }

\seealso{ See also \url{http://www.json.org/} }

\examples{
json( 1:10 )
## [1,2,3,4,5,6,7,8,9,10]

json( c("apple", "line\nwith break") )
## ["apple","line\nwith break"]

json( c(a = 1, b = 5, c = 10) )
## {"a":1,"b":5,"c":10} 

json( data.frame(a = 1:5, b = 6:10) )
## {"a":[1,2,3,4,5],"b":[6,7,8,9,10]}

json( list(a = 5, b = 1:10, c = data.frame(a = 1:2, b = 3:4)) )
## {"a":5,"b":[1,2,3,4,5,6,7,8,9,10],"c":{"a":[1,2],"b":[3,4]}}

}

\keyword{ math }
