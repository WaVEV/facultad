module Main where

import Control.Monad
import System.Environment
import System.Exit

import Eval
import Lang (Comm)
import Lexer hiding (In)
import Parser --hiding (lex)
import State (initState)


-- | La interpretación de la descripción de comportamientos de los
-- programas es interpretada como una computación en la mónada de IO.
-- Le agregue que muestre el estado en que falla
run :: Omega -> IO ()
run (Term st) = return ()
run (Abort st) =  putStrLn ("El programa terminó con una falla."++(show st))
run (Out n w) = putStrLn (show n) >> run w
run (In v g) = putStr (v ++ " : ") >> getLine >>= run . g . read

-- | El programa que pega todo: lee un archivo (y falla sin remedio si
-- no existe), les pasamos el contenido al lexer y, en caso de éxito,
-- la lista de tokens al parser. Si se pudo parsear un comando,
-- entonces lo evaluamos en el estado inicial y finalmente
-- interpretamos la descripción de su comportamiento.
main :: IO ()
main = do { args <- getArgs ;
            when (null args) errNoArg ;
            source <- readFile (head args) ;
            ts <- lexerM source ;
            comm <- parse ts ;
            run $ eval comm initState
          }

-- | Promovemos el lexing a una acción monádica, saliendo del programa
-- si hay un error.
lexerM :: String -> IO [Token]
lexerM = either (errLexer . show) return . lexer

-- | Versión monádica del parser.
parse :: [Token] -> IO Comm
parse = either (errParser . show) return . parser

-- | Fallo si no se pasa un nombre de archivo.
errNoArg :: IO a
errNoArg = putStrLn "Debe pasar un nombre de archivo como argumento." >>
           exitFailure 

-- | Fallo del lexer.
errLexer :: String -> IO a
errLexer err = putStrLn ("Error en el lexer: " ++ err) >> exitFailure


-- | Fallo del parser.
errParser :: String -> IO a
errParser err = putStrLn ("Error en el parser: " ++ err) >> exitFailure
