//! Lexer for oCAS expression syntax.
//!
//! The lexer recognises a minimal CAS language: integers, identifiers,
//! arithmetic operators `+ - * / ^`, parentheses, and commas.

use logos::Logos;

/// A token in the oCAS expression language.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Logos)]
#[logos(skip r"[ \t\n\r]+")]
#[logos(error = LexError)]
pub enum Token<'a> {
    /// An integer literal, e.g. `42` or `-7`.
    #[regex(r"-?[0-9]+", |lex| lex.slice().parse::<i64>())]
    Integer(i64),

    /// An identifier (variable or function name), e.g. `x` or `sin`.
    #[regex(r"[a-zA-Z_][a-zA-Z0-9_]*", |lex| lex.slice())]
    Ident(&'a str),

    /// `+`
    #[token("+")]
    Plus,

    /// `-`
    #[token("-")]
    Minus,

    /// `*`
    #[token("*")]
    Star,

    /// `/`
    #[token("/")]
    Slash,

    /// `^`
    #[token("^")]
    Caret,

    /// `(`
    #[token("(")]
    LParen,

    /// `)`
    #[token(")")]
    RParen,

    /// `,`
    #[token(",")]
    Comma,

    /// End-of-file sentinel added by [`lex`].
    Eof,
}

/// Lex an input string into a vector of tokens.
///
/// # Errors
///
/// Returns the first lexing error encountered.
pub fn lex(input: &str) -> Result<Vec<Token<'_>>, LexError> {
    Token::lexer(input)
        .collect::<Result<Vec<_>, _>>()
        .map(|mut tokens| {
            tokens.push(Token::Eof);
            tokens
        })
}

/// Lexing error produced when input does not match any known token.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct LexError;

impl std::fmt::Display for LexError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "invalid token")
    }
}

impl std::error::Error for LexError {}

impl From<std::num::ParseIntError> for LexError {
    fn from(_: std::num::ParseIntError) -> Self {
        LexError
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lex_integer() {
        let tokens = lex("42").unwrap();
        assert_eq!(tokens, vec![Token::Integer(42), Token::Eof]);
    }

    #[test]
    fn lex_negative_integer() {
        let tokens = lex("-7").unwrap();
        assert_eq!(tokens, vec![Token::Integer(-7), Token::Eof]);
    }

    #[test]
    fn lex_identifier() {
        let tokens = lex("x").unwrap();
        assert_eq!(tokens, vec![Token::Ident("x"), Token::Eof]);
    }

    #[test]
    fn lex_operators() {
        let tokens = lex("+ - * / ^").unwrap();
        assert_eq!(
            tokens,
            vec![
                Token::Plus,
                Token::Minus,
                Token::Star,
                Token::Slash,
                Token::Caret,
                Token::Eof
            ]
        );
    }

    #[test]
    fn lex_punctuation() {
        let tokens = lex("(, )").unwrap();
        assert_eq!(
            tokens,
            vec![Token::LParen, Token::Comma, Token::RParen, Token::Eof]
        );
    }

    #[test]
    fn lex_expression() {
        let tokens = lex("x + 2*y^3").unwrap();
        assert_eq!(
            tokens,
            vec![
                Token::Ident("x"),
                Token::Plus,
                Token::Integer(2),
                Token::Star,
                Token::Ident("y"),
                Token::Caret,
                Token::Integer(3),
                Token::Eof
            ]
        );
    }

    #[test]
    fn lex_skips_whitespace() {
        let tokens = lex("  x   \n\t+  1  ").unwrap();
        assert_eq!(
            tokens,
            vec![
                Token::Ident("x"),
                Token::Plus,
                Token::Integer(1),
                Token::Eof
            ]
        );
    }

    #[test]
    fn lex_rejects_invalid_input() {
        assert!(lex("@").is_err());
    }
}
