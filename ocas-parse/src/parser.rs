//! Recursive-descent parser for oCAS expressions.
//!
//! Parses the token stream produced by [`crate::lexer`] into an
//! [`ocas_atom::Atom`] expression tree allocated in the provided arena.

use ocas_atom::{Atom, AtomArena};
use thiserror::Error;

use crate::lexer::{LexError, Token, lex};

/// Errors that can occur while parsing an expression.
#[derive(Debug, Error, PartialEq)]
pub enum ParseError {
    /// The input could not be lexed.
    #[error("lex error")]
    Lex(#[from] LexError),
    /// Unexpected end of input.
    #[error("unexpected end of input")]
    UnexpectedEof,
    /// An unexpected token was encountered.
    #[error("unexpected token")]
    UnexpectedToken,
}

/// Parse an expression string into an [`Atom`].
///
/// # Errors
///
/// Returns a [`ParseError`] if the input cannot be lexed or parsed.
///
/// # Example
///
/// ```
/// use ocas_atom::AtomArena;
/// use ocas_core::arena::Arena;
/// use ocas_parse::parse;
///
/// let arena = Arena::new();
/// let ctx = AtomArena::new(&arena);
/// let expr = parse(&ctx, "x^2 + 2*x + 1").unwrap();
/// assert_eq!(expr.to_string(), "((x^2) + (2*x)) + 1");
/// ```
pub fn parse<'a>(ctx: &'a AtomArena<'a>, input: &'a str) -> Result<Atom<'a>, ParseError> {
    let tokens = lex(input)?;
    let mut parser = Parser::new(ctx, &tokens);
    parser.parse()
}

struct Parser<'a, 'tokens> {
    ctx: &'a AtomArena<'a>,
    tokens: &'tokens [Token<'tokens>],
    pos: usize,
}

impl<'a, 'tokens> Parser<'a, 'tokens> {
    fn new(ctx: &'a AtomArena<'a>, tokens: &'tokens [Token<'tokens>]) -> Self {
        Self {
            ctx,
            tokens,
            pos: 0,
        }
    }

    fn parse(&mut self) -> Result<Atom<'a>, ParseError> {
        let expr = self.expr()?;
        self.expect(Token::Eof)?;
        Ok(expr)
    }

    fn current(&self) -> Option<Token<'_>> {
        self.tokens.get(self.pos).copied()
    }

    fn advance(&mut self) -> Token<'_> {
        let token = self.tokens[self.pos];
        self.pos += 1;
        token
    }

    fn expect(&mut self, expected: Token) -> Result<(), ParseError> {
        match self.current() {
            Some(token) if token == expected => {
                self.advance();
                Ok(())
            }
            Some(_) => Err(ParseError::UnexpectedToken),
            None => Err(ParseError::UnexpectedEof),
        }
    }

    // expr -> term ((+|-) term)*
    fn expr(&mut self) -> Result<Atom<'a>, ParseError> {
        let mut left = self.term()?;
        while let Some(token) = self.current() {
            match token {
                Token::Plus => {
                    self.advance();
                    let right = self.term()?;
                    left = self.ctx.add(&[left, right]);
                }
                Token::Minus => {
                    self.advance();
                    let right = self.term()?;
                    let neg_right = self.ctx.mul(&[self.ctx.num(-1), right]);
                    left = self.ctx.add(&[left, neg_right]);
                }
                _ => break,
            }
        }
        Ok(left)
    }

    // term -> factor ((*|/) factor)*
    fn term(&mut self) -> Result<Atom<'a>, ParseError> {
        let mut left = self.factor()?;
        while let Some(token) = self.current() {
            match token {
                Token::Star => {
                    self.advance();
                    let right = self.factor()?;
                    left = self.ctx.mul(&[left, right]);
                }
                Token::Slash => {
                    self.advance();
                    let right = self.factor()?;
                    let neg_one = self.ctx.num(-1);
                    let inv_right = self.ctx.pow(right, neg_one);
                    left = self.ctx.mul(&[left, inv_right]);
                }
                _ => break,
            }
        }
        Ok(left)
    }

    // factor -> primary (^ factor)? with a leading minus binding looser
    // than exponentiation: -x^2 = -(x^2).
    fn factor(&mut self) -> Result<Atom<'a>, ParseError> {
        if self.current() == Some(Token::Minus) {
            self.advance();
            let operand = self.factor()?;
            return Ok(self.ctx.mul(&[self.ctx.num(-1), operand]));
        }
        let base = self.primary()?;
        if self.current() == Some(Token::Caret) {
            self.advance();
            let exp = self.factor()?;
            Ok(self.ctx.pow(base, exp))
        } else {
            Ok(base)
        }
    }

    // primary -> number | ident | ident ( arg_list ) | ( expr )
    fn primary(&mut self) -> Result<Atom<'a>, ParseError> {
        match self.current() {
            Some(Token::Integer(n)) => {
                self.advance();
                Ok(self.ctx.num(n))
            }
            Some(Token::Ident(name)) => {
                let name = name.to_owned();
                self.advance();
                if self.current() == Some(Token::LParen) {
                    self.advance();
                    let args = if self.current() == Some(Token::RParen) {
                        Vec::new()
                    } else {
                        self.arg_list()?
                    };
                    self.expect(Token::RParen)?;
                    Ok(self.ctx.fun(&name, &args))
                } else {
                    Ok(self.ctx.var(&name))
                }
            }
            Some(Token::LParen) => {
                self.advance();
                let inner = self.expr()?;
                self.expect(Token::RParen)?;
                Ok(inner)
            }
            Some(_) => Err(ParseError::UnexpectedToken),
            None => Err(ParseError::UnexpectedEof),
        }
    }

    // arg_list -> expr (',' expr)*
    fn arg_list(&mut self) -> Result<Vec<Atom<'a>>, ParseError> {
        let mut args = vec![self.expr()?];
        while self.current() == Some(Token::Comma) {
            self.advance();
            args.push(self.expr()?);
        }
        Ok(args)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ocas_core::arena::Arena;

    #[test]
    fn parse_number() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let atom = parse(&ctx, "42").unwrap();
        assert_eq!(atom.to_string(), "42");
    }

    #[test]
    fn parse_variable() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let atom = parse(&ctx, "x").unwrap();
        assert_eq!(atom.to_string(), "x");
    }

    #[test]
    fn parse_addition() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let atom = parse(&ctx, "x + y").unwrap();
        assert_eq!(atom.to_string(), "x + y");
    }

    #[test]
    fn parse_multiplication() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let atom = parse(&ctx, "x * y").unwrap();
        assert_eq!(atom.to_string(), "x*y");
    }

    #[test]
    fn parse_function_call() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let atom = parse(&ctx, "sin(x)").unwrap();
        assert_eq!(atom.to_string(), "sin(x)");
    }

    #[test]
    fn parse_function_call_multiple_args() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let atom = parse(&ctx, "f(x, y, 2)").unwrap();
        assert_eq!(atom.to_string(), "f(x, y, 2)");
    }

    #[test]
    fn parse_function_call_in_expression() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let atom = parse(&ctx, "sin(x) + cos(x)").unwrap();
        assert_eq!(atom.to_string(), "(sin(x)) + (cos(x))");
    }

    #[test]
    fn parse_operator_precedence() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let atom = parse(&ctx, "x + 2 * y").unwrap();
        assert_eq!(atom.to_string(), "x + (2*y)");
    }

    #[test]
    fn parse_right_associative_power() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let atom = parse(&ctx, "2 ^ 3 ^ 2").unwrap();
        assert_eq!(atom.to_string(), "2^(3^2)");
    }

    #[test]
    fn parse_parentheses() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let atom = parse(&ctx, "(x + y) * z").unwrap();
        assert_eq!(atom.to_string(), "(x + y)*z");
    }

    #[test]
    fn parse_unary_minus() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let atom = parse(&ctx, "-x").unwrap();
        assert_eq!(atom.to_string(), "-1*x");
    }

    #[test]
    fn parse_subtraction_normalized() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let atom = parse(&ctx, "x - y").unwrap();
        assert_eq!(atom.to_string(), "x + (-1*y)");
    }

    #[test]
    fn parse_division_normalized() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let atom = parse(&ctx, "x / y").unwrap();
        assert_eq!(atom.to_string(), "x*(y^-1)");
    }

    #[test]
    fn parse_polynomial_like() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let atom = parse(&ctx, "x^2 + 2*x + 1").unwrap();
        assert_eq!(atom.to_string(), "((x^2) + (2*x)) + 1");
    }

    #[test]
    fn parse_rejects_invalid_input() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        assert!(parse(&ctx, "x +").is_err());
    }
}
