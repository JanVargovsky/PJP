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

    public override string ToString()
    {
        return Type.ToString() + (HasValue ? $":{Value}" : "");
    }
}

class LexicalAnalyser : IDisposable
{
    private readonly TextReader input;
    private readonly Token[] operators, dividers, keywords;
    private readonly Token[] allStaticTokens;

    //private readonly char[] operators = { '+', '-', '*', '/' };
    //private readonly char[] dividers = { '(', ')', ';' };
    //private readonly string[] keywords = { "div", "mod" };
    //private readonly string[] allStaticTokens;

    public LexicalAnalyser(TextReader input)
    {
        this.input = input;
        Initialize(out operators, out dividers, out keywords);
        allStaticTokens = operators.Concat(dividers).Concat(keywords).ToArray();
        //allStaticTokens = operators.Concat(dividers.Select(t => t.ToString()).Concat(keywords.Select(t => t.ToString()).ToArray();
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

        for (int i = 0; i < probablyTokens.Length; i++)
        {
            string probablyToken = probablyTokens[i];
            // search in static tokens
            var token = allStaticTokens.FirstOrDefault(t => t.Value == probablyToken);
            if (token != null)
                yield return token;

            // number
            if (probablyToken.All(t => Char.IsDigit(t)))
                yield return new Token(TokenType.Number, probablyToken);

            StringBuilder sb = new StringBuilder();
            for (int j = 0; j < probablyToken.Length; j++)
            {
                if (probablyToken[j] == '\'')
                {
                    yield return new Token(TokenType.Comment, probablyToken.Substring(j, probablyToken.Length - j) + string.Join("", probablyTokens.Skip(i)));
                }
            }
        }
    }

    /// <summary>
    /// returns all tokens from the input
    /// </summary>
    /// <returns></returns>
    public IEnumerable<Token> GetTokens()
    {
        using (input)
        {
            string line;
            while (!string.IsNullOrEmpty((line = input.ReadLine())))
                foreach (var token in GetTokens(line))
                    yield return token;
        }
    }

    public void Dispose()
    {
        input?.Dispose();
    }
}


class Program
{
    private static TextReader GetFakeInputStream()
    {
        var ms = new MemoryStream();

        using (var sw = new StreamWriter(ms, Encoding.Default, 1024, true) { AutoFlush = true })
        {
            sw.WriteLine("    - '2 + (245 div 3);  ' poznamka");
            sw.WriteLine("    -2 + (245 div 3);  ' poznamka");
            sw.Flush();
        }

        ms.Position = 0;

        return new StreamReader(ms);
    }

    static void Main(string[] args)
    {
        using (var input = GetFakeInputStream())
        using (var lexer = new LexicalAnalyser(input))
            foreach (var item in lexer.GetTokens())
                Console.WriteLine(item);
    }
}
