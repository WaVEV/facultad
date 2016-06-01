module Examples where

import Lang

-- Ejemplos
-- x + 1
equismasuno :: IntExpr
equismasuno = v "x" .+ c 1

-- x + 1 != x + 1
assex = equismasuno .!= equismasuno


-- ejemplo de comando
commex = Newvar "x" {- := -} (c 3 .+ c 1) {- in -}
           (If (v "x" .< equismasuno)
{- then -} ("x" .:= equismasuno .# ("x" .:= equismasuno))
{- else -} ("y" .:= equismasuno))
         
commex2 = Morite

commex3 = Seq ( 
			Assign "x" (Const 1) 
		) 
		(
			Seq 
				(
					While (
						IntAss Eq (Var "y") (Const 0)
					) (
						Seq (
							If (IntAss Eq (IBin Mod (Var "x") (Const 2)) (Const 0)) ( 
								Assign "y" (IBin Plus (Var "y") (Const 1))) Skip) 
							(Seq ( 
								Assign "x" (IBin Plus (Var "x") (Const 1))
							) (Toma (IBin Minus (Var "x") (Const 0)))))) 
			Morite)

commex4 = Seq ( 
			Assign "x" (Const 1) 
		) 
		(
			While (
				IntAss Eq (Var "y") (Const 0)
			) (
				Seq (
					If (IntAss Eq (IBin Mod (Var "x") (Const 2)) (Const 0)) ( 
						Assign "y" (IBin Plus (Var "y") (Const 1))) Morite) 
					(Seq ( 
						Assign "x" (IBin Plus (Var "x") (Const 1))
					) (Toma (IBin Minus (Var "x") (Const 0)))))) 
