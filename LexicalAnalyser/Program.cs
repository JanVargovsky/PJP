using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

enum TokenType
{
    Identifier,
    Number,
    Operator,
    Divider,
    Keyword,
    Comment
}

class Token
{
    public TokenType Type { get; }
    public string Value { get; }
    public bool HasValue => !string.IsNullOrWhiteSpace(Value);

    public Token(TokenType type, string value = null)
    {
        Type = type;
        Value = value == null ? "" : value;
    }

    public override int GetHashCode()
    {
        return Type.GetHashCode() + Value.GetHashCode();
    }

    public override bool Equals(object obj)
    {
        var token = (Token)obj;
        if (token == null)
            return false;

        return Type == token.Type && Value == token.Value;
    }

    private static TokenType[] excludes = { TokenType.Divider, TokenType.Identifier, TokenType.Keyword, TokenType.Operator };

    public override string ToString()
    {
        if (excludes.Contains(Type))
            return Value;
        else
            return Type.ToString() + (HasValue ? $":{Value}" : "");
    }
}

class LexicalAnalyser
{
    //private readonly Token[] operators, dividers, keywords;
    //private readonly Token[] allStaticTokens;

    private readonly char[] operators = { '+', '-', '*', '/' };
    private readonly char[] dividers = { '(', ')', ';' };
    private readonly string[] keywords = { "div", "mod" };
    private readonly string[] allStaticTokens;

    public LexicalAnalyser()
    {
        Initialize(out operators, out dividers, out keywords);
        //allStaticTokens = operators.Concat(dividers).Concat(keywords).ToArray();
        allStaticTokens = keywords
            .Concat(dividers.Select(t => t.ToString()))
            .Concat(operators.Select(t => t.ToString()))
            .ToArray();
    }

    private void Initialize(out char[] operators, out char[] dividers, out string[] keywords)
    {
        operators = new char[] { '+', '-', '*', '/' };
        dividers = new char[] { '(', ')', ';' };
        keywords = new string[] { "div", "mod" };
    }

    private void Initialize(out Token[] operators, out Token[] dividers, out Token[] keywords)
    {
        operators = new Token[]
        {
            new Token(TokenType.Operator, "+"),
            new Token(TokenType.Operator, "-"),
            new Token(TokenType.Operator, "*"),
            new Token(TokenType.Operator, "/"),
        };

        dividers = new Token[]
        {
            new Token(TokenType.Divider, "("),
            new Token(TokenType.Divider, ")"),
            new Token(TokenType.Divider, ";"),
        };

        keywords = new Token[]
        {
            new Token(TokenType.Keyword, "div"),
            new Token(TokenType.Keyword, "mod"),
        };
    }

    public IEnumerable<Token> GetTokens(string line)
    {
        var probablyTokens = line.Split(new char[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);

        bool comment = false;
        for (int i = 0; i < probablyTokens.Length && !comment; i++)
        {
            string probablyToken = probablyTokens[i];

            if (keywords.Contains(probablyToken))
            {
                yield return new Token(TokenType.Keyword, probablyToken);
                continue;
            }

            // no spaces, so have to parse it manually
            StringBuilder identifier = new StringBuilder();
            for (int j = 0; j < probablyToken.Length && !comment; j++)
            {
                Token token = null;

                if (operators.Contains(probablyToken[j]))
                    token = new Token(TokenType.Operator, probablyToken[j].ToString());
                else if (dividers.Contains(probablyToken[j]))
                    token = new Token(TokenType.Divider, probablyToken[j].ToString());
                else if (probablyToken[j] == '\'')
                {
                    // get rest of the probablyTokens, because its comment till end of the line
                    // take comment without ' + next "tokens"
                    // TODO BUG: Take comment properly
                    token = new Token(TokenType.Comment, probablyToken.Substring(j + 1) + string.Join(" ", probablyTokens.Skip(i + 1)));
                    comment = true;
                }
                // number?
                else if (char.IsDigit(probablyToken[j]))
                {
                    int start = j;
                    while (++j < probablyToken.Length && char.IsDigit(probablyToken[j])) ;

                    token = new Token(TokenType.Number, probablyToken.Substring(start, j - start));
                    j--;
                }

                if (token != null)
                {
                    if (identifier.Length != 0)
                    {
                        yield return new Token(TokenType.Identifier, identifier.ToString());
                        identifier.Clear();
                    }
                    yield return token;
                }
                else
                    identifier.Append(probablyToken[j]);
            }

            if (identifier.Length != 0)
            {
                yield return new Token(TokenType.Identifier, identifier.ToString());
                identifier.Clear();
            }
        }
    }
}


class Program
{
    private static StreamReader GetFakeInputStream()
    {
        var ms = new MemoryStream();

        using (var sw = new StreamWriter(ms, Encoding.Default, 1024, true) { AutoFlush = true })
        {
            //sw.WriteLine("    - '2 + (245 div 3);  ' poznamka");
            //sw.WriteLine("    -2 + (245 div 3);  ' poznamka");
            sw.WriteLine(@"    - '2 + (245 div 3);  ' poznamka
- 2 + 3*(245 div 17); ' a poznamka
alfa-18*(-beta);
b*b - 4 * a * c;
' posledni poznamka
-a*b/c;");
        }

        ms.Position = 0;

        return new StreamReader(ms);
    }

    static void Main(string[] args)
    {
        var lexer = new LexicalAnalyser();

        using (var input = GetFakeInputStream())
        {
            string line;
            while (!input.EndOfStream)
            {
                line = input.ReadLine();
                Console.WriteLine($"Line: {line}");
                foreach (var item in lexer.GetTokens(line))
                    Console.WriteLine(item);
            }
        }
    }
}
