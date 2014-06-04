module Lexer where

import Data.Char(digitToInt)
import Text.Parsec

type Lexer = Parsec String ()

-- | Tokens para nuestro parser
data Token = Kw Keyword
           | Num Int
           | Op Operator
           | Str String
           | Sym Symbol
             deriving Show

-- | Símbolos
data Symbol = OpenParen  -- ^ (
            | CloseParen -- ^ )
            | OpenCurly  -- ^ {
            | CloseCurly -- ^ }
			| Admiration
			| Interrogation
              deriving (Eq,Show)

-- | Operadores
data Operator = TokLess 
              | TokEqual 
              | TokGreater
              | TokAssign
              | TokLtEq     
              | TokNot 
              | TokMinus
              | TokGtEq     
              | TokPlus 
              | TokTimes 
              | TokDiv 
              | TokMod
              | TokNotEq    
              | TokOr        
              | TokAnd             
              | TokImpl  
              | TokIff
              | TokSemicolon
      deriving (Eq,Show)
   -- | Palabras reservadas
data Keyword = KIf | Then | Else | Fi 
             | KNewVar | In 
             | KSkip | KWhile | Do | Od 
			 | KTrue | KFalse |KFail| Catchin
                deriving (Eq,Show)

-- | Lista de palabras reservadas y su representación
keywords :: [(Keyword,String)]
keywords = [ (KIf, "if")
           , (Then, "then")
           , (Else, "else")
           , (Fi, "fi")
           , (KNewVar, "newvar")
           , (In, "in")
           , (KSkip, "skip")
           , (KWhile, "while")
           , (Do, "do")
           , (Od, "od")
           , (KTrue, "true")
           , (KFalse, "false")
           , (KFail , "fail")
           , (Catchin, "catchin")
           ]

-- | Una string es una palabra reservada
reserved :: String -> Bool
reserved w = w `elem` map snd keywords

-- | Lista de operadores 
operators :: [(Operator,String)]
operators = [ (TokSemicolon, ";")
            , (TokAssign, ":=")
            , (TokIff, "<=>")
            , (TokImpl, "=>")            
            , (TokLtEq, "<=")
            , (TokGtEq, ">=")
            , (TokNotEq, "/=")                        
            , (TokLess, "<")
            , (TokEqual, "=")
            , (TokGreater, ">")
            , (TokNot, "¬")
            , (TokOr, "\\/")
            , (TokAnd, "/\\")
            , (TokMod, "%")
            , (TokMinus, "-")
            , (TokPlus, "+")
            , (TokTimes, "*")
            , (TokDiv, "/")
            ]

-- | Símbolos que no son usados para operadores. 
symbols :: [(Symbol,Char)]
symbols = [(OpenParen,'('), (CloseParen,')'), (OpenCurly,'{'), (CloseCurly,'}'), (Admiration , '!'), (Interrogation, '?')]

-- | Nuestro lexer es un parser que consume strings y devuelve
-- o bien un error o bien una lista de tokens
lexer :: String -> Either ParseError [Token]
lexer = runParser lexemes () "input"

-- | Nuestro lexer falla si no puede consumir toda la entrada.
lexemes = do _ <- spaces
             ls <- choice [ lexChoice kw keywords
                          , lexChoice op operators
                          , lexChoice sym symbols
                          , iden
                          , natural
                          ] `endBy` spaces
             _ <- eof             
             return ls

-- | Reconocimiento de un identificador
iden :: Lexer Token
iden = do var <- name
          if reserved var
          then fail "Un identificador no puede ser una palabra reservada"
          else return $ Str var
    where name :: Lexer String
          name = do h <- lower <|> char '_'
                    t <- many (alphaNum <|> char '_')
                    return (h:t)

natural :: Lexer Token
natural = do ds <- many1 digit
             return $ Num (toInt 0 ds)
    where toInt :: Int -> String -> Int
          toInt x [] = x
          toInt x (d:ds) = toInt (x*10+digitToInt d) ds

-- | @lex@ es una función que dado un par @(s,c)@ donde construye
-- un lexer para @s@ a partir de la representación @c@; @lexChoice@
-- tiene éxito con el primero de esos parser 
lexChoice :: ((a,b) -> Lexer Token) -> [(a,b)] -> Lexer Token
lexChoice lex cs = choice $ map (try . lex) cs

-- | Generación de lexers para palabras reservadas
kw :: (Keyword,String) -> Lexer Token
kw (kw,rep) = do w <- many1 letter
                 if rep == w
                 then return $ Kw kw
                 else fail "No es el keyword esperado"

-- | Generación de lexers para operadores
op :: (Operator,String) -> Lexer Token
op (o,rep) = string rep >> return (Op o)

-- | Generación de lexers para símbolos
sym :: (Symbol,Char) -> Lexer Token
sym (s,c) = do _ <- char c
               return $ Sym s 

