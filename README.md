# Ciclos econômicos e emissão de CO2 no Brasil: uma análise dinâmica para políticas ambientais ótimas



## Resumo: 

Este artigo estuda como as políticas ambientais devem responder de forma ótima aos ciclos econômicos no Brasil. Para isso utilizamos um modelo de ciclos reais de negócios (RBC) com externalidades de poluição. Os resultados indicam que o custo de mitigação de emissões de carbono é mais baixo que o da poluição no Brasil, justificando a adoção de políticas restritivas. Além disso, diante de choques de produtividade, o comportamento ótimo é uma política de mitigação pró-ciclica, onde o imposto e o teto de emissões devem aumentar em períodos de expansão da economia e diminuir em períodos de crise. Também definimos uma regra de decisão para a dinâmica da política de mitigação.

**Palavras-chave:** Políticas Ambientais; Emissão de Carbono; Ciclos Reais de Negócios; Brasil.

***

*Simulações e análises com uso do Matlab.*

***

O Brasil é signatário de diversos acordos multilaterais no esforço global de redução dos gases de efeito estufa (GEE), mesmo não tendo metas obrigatórias de mitigação de dióxido de carbono (CO2). Porém, as políticas ambientais adotadas até o momento não incluem o controle formal das emissões, tais como impostos ou comércio de quotas. A preocupação com os efeitos das mudanças climáticas tem levado alguns setores do governo a discutir a adoção de mecanismos tributários para uma economia de baixo carbono, sinalizando uma mudança de atitude. Assim, analisar como tais políticas ambientais devem responder de forma ótima aos ciclos econômicos no Brasil é essencial para auxiliar a adoção e formulação dessas políticas.

Questões ambientais como as emissões antropogênicas de GEE estão associadas a flutuações econômicas e choques de produtividade (BAUMOL e OATES, 1988). Recentemente, alguns estudos passaram a empregar modelos estocásticos dinâmicos de equilíbrio geral através da estrutura básica do modelo de Ciclos Reais de Negócios (RBC), adicionando a poluição como uma externalidade em seus modelos (FISCHER e HEUTEL, 2013). Entretanto, nenhum trabalho empregando tal metodologia foi feito para o Brasil. Assim, o objetivo deste artigo é estudar políticas ambientais formais de controle de emissão de CO2 que respondam de forma dinâmica aos ciclos econômicos no Brasil. Mais especificamente, avaliamos se (i) as reduções das externalidades geradas pelas emissões de CO2 compensam o custo desta mitigação; (ii) como as políticas de imposto e quota de emissões devem responder a choques de produtividade na economia; e (iii) qual a regra de decisão ótima que deve ser adotada para as políticas de mitigação. Para isso utilizamos um modelo RBC com externalidades de poluição proposto por Heutel (2012), calibrado com dados da economia brasileira durante o período de 1980 a 2010.

Os principais resultados apontam para três fatos. Primeiro, podemos concluir que é vantajoso adotar políticas restritivas de emissões, isto é, o custo da mitigação é inferior ao impacto da poluição sobre a produção econômica, de modo que, no equilíbrio, há uma taxa de imposto positiva que impõe um nível de mitigação maior que zero. Em segundo lugar, o comportamento ótimo do governo implica em uma política de mitigação pró-cíclica. Durante os ciclos de expansão, o nível ótimo de emissões deve crescer, mas não na mesma magnitude que cresceria sem uma política dinâmica de impostos ou quotas. Por esse motivo, tanto o imposto quanto o teto de emissões devem aumentar em períodos de expansão da economia e diminuir em períodos de crise. Por fim, um choque positivo de produtividade no modelo utilizado, por um lado aumenta o custo de mitigação, e por outro lado aumenta a demanda por mitigação. A trajetória ótima de mitigação é obtida através do aumento do imposto sobre as emissões até o ponto em que o custo de mitigação acaba anulando o aumento da demanda por mitigação. Assim, o imposto deve aumentar até o ponto em que a produtividade marginal do capital começa a cair em relação ao seu nível de estado estacionário. A partir desse momento, o comportamento ótimo do governo é o de reduzir o imposto para continuar equilibrando os custos de mitigação e emissão.

O artigo é o primeiro a estudar a adoção de políticas ambientais dinâmicas para o Brasil através de um modelo de RBC com externalidade de poluição. Os resultados podem subsidiar as discussões sobre políticas de mitigação de externalidades no Brasil, além de contribuir com a literatura de crescimento econômico.

FIGURA 2: Respostas ao impulso de +1% em $a_(t=1)$ – Variáveis econômicas
![image](https://github.com/raguirreleal/rbc_macro_ambiental--Matlab/assets/144735714/0f42e7c0-e5b7-49a0-85c4-fb4bdf7f4a2a)

FIGURA 3: Respostas ao impulso de +1% em $a_(t=1)$ – Variáveis ambientais e de política
![image](https://github.com/raguirreleal/rbc_macro_ambiental--Matlab/assets/144735714/f8fb6301-7186-4a27-b972-e26ec5b03969)

FIGURA 4: Simulação de choques aleatórios em $a_t$
![image](https://github.com/raguirreleal/rbc_macro_ambiental--Matlab/assets/144735714/e6a66348-1185-4385-8782-e191a8e6ccf4)

Os resultados nos indicam que o equilíbrio ótimo é uma taxa de imposto positiva que tenha um caráter dinâmico, aumentando em períodos de expansão econômica, e diminuindo em períodos de crise. Entretanto, o aumento da taxa de imposto em períodos de expansão não deve ser tão grande a ponto de reduzir o nível de emissões. É natural que se permita um nível de emissões maior após um choque positivo de produtividade. Isso faz com que, equivalentemente, a política de teto deva ser flexibilizada após este choque. Assim, ambas as políticas tem um caráter pró-cíclico. 

Se por um lado um choque positivo de produtividade leva a uma maior demanda por ambiente limpo, por outro lado, este choque aumenta a remuneração do capital, que é o custo de oportunidade de se investir em políticas de mitigação. Para o caso do Brasil, vimos que o primeiro efeito domina o segundo, de modo que devemos aumentar a taxa de imposto após o choque. Porém, no momento em que a taxa de remuneração do capital começar a variar negativamente, isso implica que a produtividade marginal do capital está caindo, de modo que a imposição da alíquota de imposto começa a penalizar a produção além do ponto ótimo, e então o imposto deve começar a cair. Esta é a regra ótima para a definição da política dinâmica.


