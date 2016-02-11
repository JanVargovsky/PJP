using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

internal enum AutomatRunResult
{
    Final, NotFinal, NotDeterministic,
}

internal class Automat
{
    public List<Transition> Transitions { get; set; }
    public string InputSymbols { get; set; }
    public List<int> FinalStates { get; set; }

    public Automat()
    {
        Transitions = new List<Transition>();
    }

    public bool Run(string input)
    {
        if (string.IsNullOrEmpty(input))
            return FinalStates.Contains(0);

        int currentState = 0, nextState;

        foreach (var symbol in input)
        {
            nextState = GetNextState(currentState, symbol);
            if (nextState == -1)

            {
                Console.WriteLine($"neni deterministicky pri stavu {currentState} a symbolu {symbol}");
                return false;
            }

            currentState = nextState;
        }

        return FinalStates.Contains(currentState);
    }

    private int GetNextState(int currentState, char currentSymbol)
    {
        var state = Transitions.Find(s => s.InputState == currentState && s.Symbol == currentSymbol);

        return state == null ? -1 : state.OutputState;
    }
}

internal class Transition
{
    public int InputState { get; set; }
    public char Symbol { get; set; }
    public int OutputState { get; set; }

    public override string ToString()
    {
        return $"{InputState} {Symbol} {OutputState}";
    }
}

internal class Program
{

    static void Main(string[] args)
    {
        Console.SetIn(new StreamReader(@"..\..\test2.txt"));

        int casesCount = int.Parse(Console.ReadLine());
        while (casesCount-- > 0)
        {
            var automat = new Automat();
            // states
            int statesCount = int.Parse(Console.ReadLine());
            automat.InputSymbols = Console.ReadLine();

            // symbols
            int finalStatesCount = int.Parse(Console.ReadLine());
            var states = Console.ReadLine();
            automat.FinalStates = states.Split(' ').Select(int.Parse).ToList();

            // change functions
            int edgesCount = int.Parse(Console.ReadLine());
            for (int i = 0; i < edgesCount; i++)
            {
                string line = Console.ReadLine();
                automat.Transitions.Add(new Transition
                {
                    InputState = line[0] - '0',
                    Symbol = line[2],
                    OutputState = line[4] - '0',
                });
            }


            int testsCount = int.Parse(Console.ReadLine());
            for (int i = 0; i < testsCount; i++)
            {
                var result = automat.Run(Console.ReadLine());
                Console.WriteLine(result ? "ANO" : "NE");
            }

            if (casesCount != 0)
                Console.WriteLine();
        }
    }
}