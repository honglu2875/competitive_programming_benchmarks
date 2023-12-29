use std::collections::HashMap;
use std::io;
use std::mem::swap;

fn main() -> io::Result<()> {
    let mut buffer = String::new();
    io::stdin().read_line(&mut buffer)?;
    match buffer.trim() {
        "ENCODE" => encode(),
        "DECODE" => decode(),
        _ => (),
    };
    Ok(())
}

fn encode() {
    let mut buffer = String::new();
    let mut map = HashMap::new();
    let mut ceo = None;
    while io::stdin().read_line(&mut buffer).is_ok() && buffer.trim() != "" {
        let mut parts = buffer.trim().split(": ");
        let boss: String = parts.next().unwrap().to_owned();
        let employees = parts
            .next()
            .unwrap()
            .split(" ")
            .map(|s| s.to_owned())
            .collect::<Vec<String>>();
        map.insert(boss.clone(), employees);
        if ceo.is_none() {
            ceo = Some(boss);
        }
        buffer.clear();
    }
    println!("{}", cypher(&map, &ceo.unwrap().as_str()));
}

fn cypher(map: &HashMap<String, Vec<String>>, boss: &str) -> String {
    println!("{boss}");
    match map.get(boss) {
        None => String::new(),
        Some(employees) => employees
            .iter()
            .map(|employee| format!("0{}1", cypher(map, employee)))
            .fold(String::new(), |acc, x| acc + &x),
    }
}

fn decode() {
    let mut buffer = String::new();
    let mut ceo = None;
    let mut names = Vec::new();
    let mut code = String::new();
    while io::stdin().read_line(&mut buffer).is_ok() && buffer.trim() != "" {
        if ceo.is_none() {
            ceo = Some(buffer.trim().to_owned());
        } else {
            names.push(buffer.trim().to_owned());
        }
        swap(&mut code, &mut buffer);
        buffer.clear();
    }
    names.pop(); // Pop code from the list of names
    let mut current = ceo.as_ref().unwrap().to_owned();
    let mut code = code.chars();
    let mut graph: HashMap<String, Vec<String>> = HashMap::new();
    let mut managers: HashMap<String, String> = HashMap::new();
    for name in names {
        while code.next() == Some('1') {
            current = managers.get(&current).unwrap().to_string();
        }
        graph.entry(current.clone()).or_default().push(name.clone());
        managers.insert(name.clone(), current.clone());
        current = name.clone();
    }
    let mut q = vec![ceo.unwrap()];
    let mut prevq = vec![];
    while !q.is_empty() {
        for name in &q {
            if let Some(employees) = graph.get(name) {
                prevq.append(&mut employees.clone());
                println!("{name}: {}", employees.join(" "));
            }
        }
        swap(&mut q, &mut prevq);
        prevq.clear();
    }
}
